"""Pedigree reconstruction from per-pair :RELATIVE labels.

V1 algorithm (deliberately simple — PRIMUS-lite):

1. For each family (a connected component on :RELATIVE edges from
   M5's families procedure), pull every per-pair label.
2. Identify parent-child edges. Each is a directed parent → child
   edge in the working DiGraph; orientation comes from
   :Sample.sex when present, otherwise from age-or-id heuristic.
3. Identify full-sib edges. Sibs share two parents. v1 writes both
   parents as 0 (unknown) when neither parent appears in the
   per-pair labels; otherwise it picks the lex-smallest parent
   from the parent_child edges.
4. Emit PED rows: (FID, IID, FatherID, MotherID, Sex, Phenotype).

Defers: half-siblings (need 2nd-degree disambiguation),
3rd-degree (cousins), grandparent ↔ avuncular orientation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import networkx as nx


SEX_MALE = 1
SEX_FEMALE = 2
SEX_UNKNOWN = 0
PHENOTYPE_MISSING = -9


def _coerce_sex(value) -> int:
    """Normalise sex to PED integer encoding (0/1/2)."""
    if value is None:
        return SEX_UNKNOWN
    if isinstance(value, int):
        return value if value in (SEX_MALE, SEX_FEMALE) else SEX_UNKNOWN
    if isinstance(value, str):
        v = value.strip().upper()
        if v in ("M", "MALE", "1"):
            return SEX_MALE
        if v in ("F", "FEMALE", "2"):
            return SEX_FEMALE
    return SEX_UNKNOWN


@dataclass
class PedigreeRow:
    family_id: str
    sample_id: str
    father_id: str        # "0" for unknown
    mother_id: str        # "0" for unknown
    sex: int              # 1 = male, 2 = female, 0 = unknown
    phenotype: int = PHENOTYPE_MISSING


@dataclass
class _SampleMeta:
    sample_id: str
    sex: int = SEX_UNKNOWN


@dataclass
class _RelativeEdge:
    sample_a: str
    sample_b: str
    relationship: str
    degree: int


@dataclass
class _Family:
    family_id: str
    samples: list[str] = field(default_factory=list)
    metadata: dict[str, _SampleMeta] = field(default_factory=dict)
    edges: list[_RelativeEdge] = field(default_factory=list)


class PedigreeReconstructor:
    """Reconstruct PED rows from a Neo4j driver + family selection.

    Caller supplies an iterable of (family_id, sample_id, sex) and
    a list of per-pair :RELATIVE labels. The class is graph-driver-
    agnostic so tests can build fixtures in memory.
    """

    def __init__(
        self,
        samples: Iterable[tuple[str, str, int | None]],
        relative_edges: Iterable[tuple[str, str, str, int]],
        source: str = "king",
    ):
        """Build the in-memory family map.

        samples         : iterable of (family_id, sample_id, sex_int_or_None)
        relative_edges  : iterable of (sample_a, sample_b, relationship, degree)
        """
        self.source = source
        self._families: dict[str, _Family] = {}
        for family_id, sid, sex in samples:
            fam = self._families.setdefault(family_id, _Family(family_id))
            fam.samples.append(sid)
            fam.metadata[sid] = _SampleMeta(
                sid, _coerce_sex(sex))
        sample_to_family: dict[str, str] = {
            sid: fam_id
            for fam_id, fam in self._families.items()
            for sid in fam.samples
        }
        for a, b, rel, deg in relative_edges:
            fa = sample_to_family.get(a)
            fb = sample_to_family.get(b)
            if fa is None or fb is None or fa != fb:
                continue
            self._families[fa].edges.append(
                _RelativeEdge(a, b, rel, int(deg)))

    def reconstruct_family(self, family_id: str) -> list[PedigreeRow]:
        fam = self._families.get(family_id)
        if fam is None:
            return []
        return self._reconstruct(fam)

    def reconstruct_all(self) -> list[PedigreeRow]:
        rows: list[PedigreeRow] = []
        for fam in self._families.values():
            rows.extend(self._reconstruct(fam))
        return rows

    def _reconstruct(self, fam: _Family) -> list[PedigreeRow]:
        # Build undirected parent_child multigraph and a separate
        # full_sibling graph.
        pc = nx.MultiGraph()
        pc.add_nodes_from(fam.samples)
        sib_groups = nx.Graph()
        sib_groups.add_nodes_from(fam.samples)
        for e in fam.edges:
            if e.relationship == "parent_child":
                pc.add_edge(e.sample_a, e.sample_b)
            elif e.relationship == "full_sibling":
                sib_groups.add_edge(e.sample_a, e.sample_b)

        parents: dict[str, list[str]] = {sid: [] for sid in fam.samples}

        # Peel-off algorithm. A node with exactly 2 :parent_child
        # edges in the current graph is unambiguously a *child* —
        # its 2 neighbours are its parents. Iteratively peel these
        # leaves; multi-generation pedigrees resolve from outermost
        # leaves inward.
        progress = True
        while progress:
            progress = False
            for node in list(pc.nodes()):
                if pc.degree(node) != 2 or parents[node]:
                    continue
                ps = sorted(set(pc.neighbors(node)))
                parents[node] = ps
                for p in ps:
                    while pc.has_edge(node, p):
                        pc.remove_edge(node, p)
                progress = True

        # Remaining deg-1 :parent_child edges: orient by sex when one
        # endpoint has known sex (the M/F endpoint is the parent).
        progress = True
        while progress:
            progress = False
            for u, v in list(pc.edges()):
                if parents[u] or parents[v]:
                    continue
                u_sex = self._has_sex(u, fam.metadata)
                v_sex = self._has_sex(v, fam.metadata)
                if u_sex and not v_sex:
                    parents[v] = [u]
                    pc.remove_edge(u, v)
                    progress = True
                    break
                if v_sex and not u_sex:
                    parents[u] = [v]
                    pc.remove_edge(u, v)
                    progress = True
                    break

        # Anything left: orient lex-smaller → parent (deterministic
        # but possibly wrong; sex metadata fixes most of these).
        for u, v in list(pc.edges()):
            if parents[u] or parents[v]:
                continue
            parent = u if u < v else v
            child = v if u < v else u
            parents[child] = [parent]
            pc.remove_edge(u, v)

        # Sib groups: full-siblings share two parents. Propagate from
        # any sibling that already has parents to the rest of its group.
        for component in nx.connected_components(sib_groups):
            if len(component) < 2:
                continue
            shared = []
            for sid in component:
                shared.extend(parents.get(sid, []))
            shared = sorted(set(shared))
            if shared:
                for sid in component:
                    parents[sid] = shared[:]

        # Resolve father / mother slots from sex metadata when possible.
        rows = []
        for sid in sorted(fam.samples):
            meta = fam.metadata[sid]
            ps = sorted(set(parents.get(sid, [])))
            father = "0"
            mother = "0"
            for p in ps:
                pmeta = fam.metadata.get(p)
                psex = pmeta.sex if pmeta else SEX_UNKNOWN
                if psex == SEX_MALE and father == "0":
                    father = p
                elif psex == SEX_FEMALE and mother == "0":
                    mother = p
                elif father == "0":
                    father = p
                elif mother == "0":
                    mother = p
            rows.append(PedigreeRow(
                family_id=fam.family_id,
                sample_id=sid,
                father_id=father,
                mother_id=mother,
                sex=meta.sex,
            ))
        return rows

    @staticmethod
    def _has_sex(sid: str, meta: dict[str, _SampleMeta]) -> bool:
        m = meta.get(sid)
        return m is not None and m.sex in (SEX_MALE, SEX_FEMALE)
