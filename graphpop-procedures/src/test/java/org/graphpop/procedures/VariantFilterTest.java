package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link VariantFilter}.
 */
class VariantFilterTest {

    // ---- isActive / inactive ----

    @Test
    void nullOptions_isInactive() {
        VariantFilter f = VariantFilter.fromOptions(null);
        assertFalse(f.isActive());
    }

    @Test
    void emptyOptions_isInactive() {
        VariantFilter f = VariantFilter.fromOptions(Map.of());
        assertFalse(f.isActive());
    }

    @Test
    void defaultValues_isInactive() {
        Map<String, Object> opts = new HashMap<>();
        opts.put("min_af", 0.0);
        opts.put("max_af", 1.0);
        VariantFilter f = VariantFilter.fromOptions(opts);
        assertFalse(f.isActive());
    }

    // ---- min_af filter ----

    @Test
    void minAf_filtersBelowThreshold() {
        VariantFilter f = VariantFilter.fromOptions(Map.of("min_af", 0.05));
        assertTrue(f.isActive());

        // af=0.03 < 0.05 → should fail
        assertFalse(f.passes(3, 100, 0.03, 3, 0, null));

        // af=0.10 > 0.05 → should pass
        assertTrue(f.passes(10, 100, 0.10, 10, 0, null));

        // af=0.05 exactly → should pass
        assertTrue(f.passes(5, 100, 0.05, 5, 0, null));
    }

    // ---- max_af filter ----

    @Test
    void maxAf_filtersAboveThreshold() {
        VariantFilter f = VariantFilter.fromOptions(Map.of("max_af", 0.95));
        assertTrue(f.isActive());

        // af=0.97 > 0.95 → should fail
        assertFalse(f.passes(97, 100, 0.97, 3, 47, null));

        // af=0.90 < 0.95 → should pass
        assertTrue(f.passes(90, 100, 0.90, 10, 40, null));
    }

    // ---- combined min_af + max_af ----

    @Test
    void combinedAfFilter() {
        Map<String, Object> opts = new HashMap<>();
        opts.put("min_af", 0.05);
        opts.put("max_af", 0.95);
        VariantFilter f = VariantFilter.fromOptions(opts);

        assertTrue(f.passes(50, 100, 0.50, 25, 12, null));  // mid-range
        assertFalse(f.passes(2, 100, 0.02, 2, 0, null));    // too rare
        assertFalse(f.passes(98, 100, 0.98, 2, 48, null));  // too common
    }

    // ---- HWE filter ----

    @Test
    void hweFilter_excludesOutOfHWE() {
        // A significant HWE departure: excess heterozygosity
        // nTotal=100, nHet=80, nHomMinor=0 → extreme excess het
        // The mid-p pvalue for this should be very small
        VariantFilter f = VariantFilter.fromOptions(Map.of("hwe_pvalue", 0.001));
        assertTrue(f.isActive());

        // Extreme excess het: 80 hets, 0 hom minor, 20 hom major
        // AC = 80, AN = 200, hetCount=80, homAlt=0
        boolean passes = f.passes(80, 200, 0.40, 80, 0, null);
        assertFalse(passes, "Extreme excess het should fail HWE filter");
    }

    @Test
    void hweFilter_passesWhenInHWE() {
        // A variant in HWE: p=0.5, expected het=0.5, n=100
        // Expected: homMajor=25, het=50, homMinor=25
        VariantFilter f = VariantFilter.fromOptions(Map.of("hwe_pvalue", 0.001));

        boolean passes = f.passes(100, 200, 0.50, 50, 25, null);
        assertTrue(passes, "Variant in HWE should pass filter");
    }

    // ---- variant_type filter ----
    // Note: variant_type requires a Node, which we don't have in unit tests.
    // We test that the filter logic parses correctly and null nodes are handled.

    @Test
    void variantType_activeWhenSet() {
        VariantFilter f = VariantFilter.fromOptions(Map.of("variant_type", "SNP"));
        assertTrue(f.isActive());

        // With null node, variant_type check is skipped (cannot access property)
        // so the variant passes on other criteria
        assertTrue(f.passes(10, 100, 0.10, 10, 0, null));
    }

    // ---- no filter = always passes ----

    @Test
    void noFilter_alwaysPasses() {
        VariantFilter f = VariantFilter.fromOptions(Map.of());

        // Even extreme values should pass when no filter active
        assertTrue(f.passes(0, 100, 0.0, 0, 0, null));
        assertTrue(f.passes(100, 100, 1.0, 0, 50, null));
    }

    // ---- parsing edge cases ----

    @Test
    void numericParsing_integerValues() {
        // Options map may contain Integer instead of Double
        Map<String, Object> opts = new HashMap<>();
        opts.put("min_af", 0);       // integer 0
        opts.put("max_af", 1);       // integer 1
        opts.put("hwe_pvalue", 0);
        VariantFilter f = VariantFilter.fromOptions(opts);
        assertFalse(f.isActive());
    }
}
