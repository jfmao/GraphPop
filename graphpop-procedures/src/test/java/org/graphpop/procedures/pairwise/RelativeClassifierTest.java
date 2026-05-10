package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class RelativeClassifierTest {

    @Test
    void identical_at_threshold() {
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.45, 0.0);
        assertEquals("identical", v.relationship);
        assertEquals(0, v.degree);
    }

    @Test
    void identical_above_threshold() {
        // MZ twins: phi very close to 0.5
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.50, 0.0);
        assertEquals("identical", v.relationship);
    }

    @Test
    void parent_child_distinguished_from_full_sibling_by_ibs0() {
        // First-degree (phi ~ 0.25)
        RelativeClassifier.Verdict pc =
                RelativeClassifier.classify(0.25, 0.001);
        assertEquals("parent_child", pc.relationship);
        assertEquals(1, pc.degree);

        RelativeClassifier.Verdict fs =
                RelativeClassifier.classify(0.25, 0.06);
        assertEquals("full_sibling", fs.relationship);
        assertEquals(1, fs.degree);
    }

    @Test
    void first_degree_label_when_ibs0_unavailable() {
        // NaN ibs0 -> can't disambiguate, return umbrella label.
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.25, Double.NaN);
        assertEquals("first_degree", v.relationship);
        assertEquals(1, v.degree);
    }

    @Test
    void second_degree() {
        // half-sib / avuncular / grandparent: phi ~ 0.125
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.125, 0.05);
        assertEquals("second_degree", v.relationship);
        assertEquals(2, v.degree);
    }

    @Test
    void third_degree() {
        // first cousin: phi ~ 0.0625
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.0625, 0.05);
        assertEquals("third_degree", v.relationship);
        assertEquals(3, v.degree);
    }

    @Test
    void unrelated_below_third_degree_threshold() {
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.02, 0.10);
        assertEquals("unrelated", v.relationship);
        assertEquals(99, v.degree);
    }

    @Test
    void boundary_cases_exact_thresholds() {
        // At each threshold, the >= cutoff branch is taken.
        // identical cutoff = 0.354
        assertEquals("identical",
                RelativeClassifier.classify(0.354, 0.0).relationship);
        // firstDegree cutoff = 0.177
        assertEquals("parent_child",
                RelativeClassifier.classify(0.177, 0.0).relationship);
        // secondDegree cutoff = 0.0884
        assertEquals("second_degree",
                RelativeClassifier.classify(0.0884, 0.0).relationship);
        // thirdDegree cutoff = 0.0442
        assertEquals("third_degree",
                RelativeClassifier.classify(0.0442, 0.0).relationship);
        // Below third-degree
        assertEquals("unrelated",
                RelativeClassifier.classify(0.0441, 0.0).relationship);
    }

    @Test
    void custom_cutoffs_take_precedence() {
        // Make every pair "identical" by lowering the bar
        RelativeClassifier.Cutoffs custom =
                new RelativeClassifier.Cutoffs(0.0, 0.0, 0.0, 0.0, 0.0);
        RelativeClassifier.Verdict v =
                RelativeClassifier.classify(0.001, 0.0, custom);
        assertEquals("identical", v.relationship);
    }

    @Test
    void custom_ibs0_threshold() {
        // Bump ibs0 threshold; now ibs0_frac=0.06 stays parent_child.
        RelativeClassifier.Cutoffs custom =
                new RelativeClassifier.Cutoffs(0.354, 0.177, 0.0884, 0.0442, 0.10);
        assertEquals("parent_child",
                RelativeClassifier.classify(0.25, 0.06, custom).relationship);
    }
}
