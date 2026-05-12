//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: bench/test_ydrop.c
//
// Validation driver for both textbook ports in lastz/src/ydrop_sane.c:
//
//   - ydrop_one_sided_align_impl_sane                      ("sane")
//   - ydrop_one_sided_align_impl_sane_double_buffered      ("sane_db")
//
// Each is compared against lastz's ydrop_one_sided_align (via the
// for_testing trampoline in gapped_extend.c) on a small corpus of cases.
//
// For each case we run the comparison in TWO directions:
//   - forward direction:  all three impls on (A, B)
//   - reverse direction:  all three impls on (rev(A), rev(B)), with
//                         lastz called in reversed=1 mode
//
// Per case this yields 4 PASS/FAIL points: {sane, sane_db} x {fwd, rev}.
//
// Exit code 0 iff every comparison matches lastz on (score, end1, end2,
// edit-script). Exit code 1 on the first FAIL with diagnostic detail.
//
//----------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "build_options.h"
#include "utilities.h"
#include "dna_utilities.h"
#include "edit_script.h"
#include "lastz.h"
#include "gapped_extend.h"
#include "ydrop_sane.h"

// ---- stubs for lastz globals normally defined in lastz.c -----------------
// These live in lastz.c (which has main() and is excluded from our link).
// output.c references them via `extern` decls. We never call any output.c
// path, so the contents don't matter; we just need symbols to satisfy the
// linker. control* currParams stays NULL — any path that dereferences it
// would crash, but nothing in our test reaches there.
char *programName            = "test_ydrop";
char *programVersionMajor    = "0";
char *programVersionMinor    = "0";
char *programVersionSubMinor = "0";
char *programRevisionDate    = "20260506";
char *svnRevisionNumber      = "0";
control *currParams          = NULL;

// gapped_extend.o references these (in gappily_extend_hsps, which we never
// call) but the linker still requires them. Stubs that suicide-loudly if
// somehow reached.
void print_align_list (FILE *f, void *a)
    { (void) f; (void) a; abort (); }
void print_align_list_segments (FILE *f, void *a, int b, char *c, char *d, int e)
    { (void) f; (void) a; (void) b; (void) c; (void) d; (void) e; abort (); }

// ---- small helpers ------------------------------------------------------

// Allocate a writable copy of a literal string. lastz takes `u8*`, not
// `const u8*`, so we hand it freshly-allocated buffers.
static u8 *dup_seq (const char *s, unspos *len_out)
    {
    unspos n = (unspos) strlen (s);
    u8 *p = (u8*) malloc_or_die ("dup_seq", n + 1);
    memcpy (p, s, n);
    p[n] = 0;
    *len_out = n;
    return p;
    }

// In-place byte reversal (NOT reverse-complement).
static u8 *reverse_seq (const u8 *src, unspos n)
    {
    u8 *r = (u8*) malloc_or_die ("reverse_seq", n + 1);
    for (unspos i = 0; i < n; i++) r[i] = src[n - 1 - i];
    r[n] = 0;
    return r;
    }

// Compare two edit scripts byte-for-byte. Both should be in the same
// canonical merged-run form because both impls call edit_script_add with
// 1-base operations, which merges consecutive same-op calls identically.
static int edit_scripts_equal (editscript *a, editscript *b)
    {
    if (a == NULL && b == NULL) return 1;
    if (a == NULL || b == NULL) return 0;
    if (a->len != b->len) return 0;
    for (u32 i = 0; i < a->len; i++)
        if (a->op[i] != b->op[i]) return 0;
    return 1;
    }

static void print_edit_script (FILE *f, const char *tag, editscript *s)
    {
    if (s == NULL) { fprintf (f, "  %s: <null>\n", tag); return; }
    fprintf (f, "  %s: len=%u  ", tag, (unsigned) s->len);
    u32 to_print = s->len < 12 ? s->len : 12;
    for (u32 i = 0; i < to_print; i++)
        {
        u32 op  = edit_op_operation (s->op[i]);
        u32 rpt = edit_op_repeat    (s->op[i]);
        const char *name = (op == editopSub) ? "S"
                        : (op == editopIns) ? "I"
                        : (op == editopDel) ? "D" : "?";
        fprintf (f, "%s%u ", name, rpt);
        }
    if (s->len > to_print) fprintf (f, "...");
    fprintf (f, "\n");
    }

// Function-pointer type matching both sane impls' signatures.
typedef score (*sane_impl_fn) (
    const u8 *, unspos,
    const u8 *, unspos,
    scorerow *,
    score, score, score,
    editscript **,
    unspos *, unspos *);

typedef struct sane_impl_descr
    {
    const char  *name;       // short tag for logs
    sane_impl_fn fn;
    } sane_impl_descr;

static sane_impl_descr SANE_IMPLS[] = {
    { "sane",    ydrop_one_sided_align_impl_sane                   },
    { "sane_db", ydrop_one_sided_align_impl_sane_double_buffered   },
};
static const int N_SANE_IMPLS =
    (int)(sizeof(SANE_IMPLS) / sizeof(SANE_IMPLS[0]));

// ---- one impl-vs-lastz comparison ---------------------------------------
//
// Compares ONE sane impl against lastz on a given (A, B, direction). Lastz
// is the reference. Returns 1 on pass, 0 on fail; prints a one-line
// PASS/FAIL record either way.
//
static int compare_one
   (const char *case_name, const char *direction,
    const sane_impl_descr *impl,
    u8 *A, unspos M, u8 *B, unspos N,
    int reversed,
    scoreset *scoring, score yDrop, tback *tb)
    {
    editscript *es_lastz = edit_script_new ();
    editscript *es_sane  = edit_script_new ();

    score  sL, sS;
    unspos e1L, e2L, e1S, e2S;

    sL = ydrop_one_sided_align_for_testing (
             A, B, M, N, scoring, yDrop, tb,
             reversed, /*trimToPeak=*/1,
             &es_lastz, &e1L, &e2L);

    sS = impl->fn (
             A, M, B, N, scoring->sub,
             scoring->gapOpen, scoring->gapExtend, yDrop,
             &es_sane, &e1S, &e2S);

    int ok = (sL == sS) && (e1L == e1S) && (e2L == e2S)
          && edit_scripts_equal (es_lastz, es_sane);

    char label[80];
    snprintf (label, sizeof(label), "%s [%s]", case_name, impl->name);

    if (ok)
        printf ("[PASS] %-40s %-7s  score=%d  end=(%lu,%lu)  scriptlen=%u\n",
                label, direction,
                (int) sL, (unsigned long) e1L, (unsigned long) e2L,
                (unsigned) es_lastz->len);
    else
        {
        printf ("[FAIL] %-40s %-7s\n", label, direction);
        printf ("  lastz:   score=%d end=(%lu,%lu)\n",
                (int) sL, (unsigned long) e1L, (unsigned long) e2L);
        printf ("  %-6s:  score=%d end=(%lu,%lu)\n",
                impl->name, (int) sS,
                (unsigned long) e1S, (unsigned long) e2S);
        print_edit_script (stdout, "lastz script", es_lastz);
        print_edit_script (stdout, "sane  script", es_sane);
        }

    free_if_valid ("test_ydrop es_lastz", es_lastz);
    free_if_valid ("test_ydrop es_sane",  es_sane);
    return ok;
    }

// Run BOTH directions x ALL sane impls on a case (via input reversal for
// the reverse test). Returns total comparisons performed and how many
// passed, via out-params.
static void run_case (const char *name, const char *Astr, const char *Bstr,
                      score yDrop, scoreset *scoring, tback *tb,
                      int *out_total, int *out_passed)
    {
    unspos M, N;
    u8 *A  = dup_seq (Astr, &M);
    u8 *B  = dup_seq (Bstr, &N);
    u8 *rA = reverse_seq (A, M);
    u8 *rB = reverse_seq (B, N);

    int total = 0, passed = 0;

    for (int i = 0; i < N_SANE_IMPLS; i++)
        {
        total++;
        if (compare_one (name, "forward",
                         &SANE_IMPLS[i],
                         A, M, B, N, /*reversed=*/0,
                         scoring, yDrop, tb))
            passed++;

        // For the reversed direction: lastz expects the ORIGINAL strings
        // A,B and reversed=1, while our sane impls get the byte-reversed
        // strings. Both should produce the same alignment shape.
        //
        // (lastz's ydrop_one_sided_align(reversed=1) is the same DP run on
        // the sequences that the caller would have prepared as rev1, rev2
        // — i.e. byte-reversed. Our sane forward-only impls on rA, rB do
        // that directly.)
        total++;
        if (compare_one (name, "reverse",
                         &SANE_IMPLS[i],
                         rA, M, rB, N, /*reversed=*/1,
                         scoring, yDrop, tb))
            passed++;
        }

    free_if_valid ("test_ydrop A",  A);
    free_if_valid ("test_ydrop B",  B);
    free_if_valid ("test_ydrop rA", rA);
    free_if_valid ("test_ydrop rB", rB);

    *out_total  = total;
    *out_passed = passed;
    }

// ---- random sequence generator (xorshift, deterministic) ----------------

static u64 rng_state = 0xdeadbeefcafebabeULL;
static void   rng_seed (u64 s) { rng_state = s ? s : 0xdeadbeefcafebabeULL; }
static u32    rng_next (void) {
    u64 x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return (u32) x;
}
static u32 rng_in (u32 lo, u32 hi) { return lo + (rng_next () % (hi - lo + 1)); }

static char rand_base (void)
    {
    static const char alphabet[] = "ACGT";
    return alphabet[rng_next () & 3];
    }

// Generate two random sequences A,B at approximately `identity` percent
// (where B is a mutation of A; insertion/deletion roughly 1% probability
// each, mismatch otherwise). The output strings are heap-allocated and
// must be free()'d by the caller.
static void gen_random_pair (unspos len_A, unsigned identity_pct,
                              char **out_A, char **out_B)
    {
    char *A = (char*) malloc_or_die ("gen A", len_A + 1);
    // Generate a length for B that is within a few bp of len_A.
    char *B = (char*) malloc_or_die ("gen B", 2 * len_A + 1);
    for (unspos i = 0; i < len_A; i++) A[i] = rand_base ();
    A[len_A] = 0;

    unspos i = 0;   // index into A
    unspos j = 0;   // index into B
    while (i < len_A)
        {
        u32 r = rng_in (0, 99);
        if (r < identity_pct - 2)        // copy as-is
            { B[j++] = A[i++]; }
        else if (r < identity_pct - 1)   // 1% chance: insertion in B
            { B[j++] = rand_base (); }
        else if (r < identity_pct)       // 1% chance: deletion in B
            { i++; }
        else                              // mismatch
            {
            char c = A[i++];
            char nc;
            do { nc = rand_base (); } while (nc == c);
            B[j++] = nc;
            }
        if (j >= 2 * len_A) break;
        }
    B[j] = 0;
    *out_A = A;
    *out_B = B;
    }

// ---- main ---------------------------------------------------------------

int main (int argc, char **argv)
    {
    (void) argc; (void) argv;

    // ---- Build a HOXD-style scoring set ----
    score gapOpen = 400, gapExtend = 30;
    scoreset *scoring = new_dna_score_set (HOXD70, HOXD70_X, HOXD70_fill,
                                            gapOpen, gapExtend);
    if (scoring == NULL) { fprintf (stderr, "scoreset alloc failed\n"); return 1; }

    // ---- Pre-allocate a traceback tape (lastz needs this for the wrapper) ----
    tback *tb = new_traceback (10 * 1024 * 1024);   // 10 MB
    if (tb == NULL) { fprintf (stderr, "traceback alloc failed\n"); return 1; }

    int total = 0, passed = 0;

    // ---- Hand-crafted cases ----
    struct hc_case { const char *name, *A, *B; score yDrop; };
    struct hc_case cases[] = {
        { "01_identity_8bp",   "ACGTACGT",   "ACGTACGT",                   910 },
        { "02_one_mismatch",   "ACGT",       "AGGT",                        910 },
        { "03_one_insertion",  "ACGT",       "ACGGT",                       910 },
        { "04_one_deletion",   "ACGGT",      "ACGT",                        910 },
        { "05_ydrop_fires",    "AAAA",       "AAAATTTTTTTTTTTTTTTTTTTT",     50 },
        { "06_alternating",    "ACATATATATATATATATATAC",
                                "ACACACACACACACACACACAC",                   910 },
        { "07_long_identity",  "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA",
                                "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA", 910 },
        { "08_short_no_match", "AAAA",       "TTTT",                         50 },
    };

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++)
        {
        int t = 0, p = 0;
        run_case (cases[i].name, cases[i].A, cases[i].B,
                  cases[i].yDrop, scoring, tb, &t, &p);
        total  += t;
        passed += p;
        }

    // ---- Fuzz loop: random pairs at varying length / identity ----
    // Default 200 cases; override via FUZZ=<N> env var. Each case produces
    // (2 directions x N_SANE_IMPLS) comparison points.
    int fuzz_runs = 200;
    {
    const char *env = getenv ("FUZZ");
    if (env != NULL && env[0] != 0)
        {
        int n = atoi (env);
        if (n >= 0) fuzz_runs = n;
        }
    }
    rng_seed (1234567);
    for (int k = 0; k < fuzz_runs; k++)
        {
        unspos len = rng_in (50, 1500);
        unsigned ident = rng_in (60, 99);
        char *A, *B;
        gen_random_pair (len, ident, &A, &B);

        char tag[64];
        snprintf (tag, sizeof(tag), "fuzz_%03d_len%u_id%u", k, (unsigned) len, ident);

        int t = 0, p = 0;
        run_case (tag, A, B, /*yDrop=*/910, scoring, tb, &t, &p);
        total  += t;
        passed += p;

        free (A);
        free (B);
        }

    // ---- Cleanup ----
    free_traceback (tb);
    free_score_set ("test_ydrop scoring", scoring);

    int failed = total - passed;
    printf ("\n=========================================\n");
    printf ("test_ydrop summary: %d / %d passed (%d failed)\n",
            passed, total, failed);
    printf ("(per case: %d sane impls x 2 directions = %d comparisons)\n",
            N_SANE_IMPLS, 2 * N_SANE_IMPLS);
    printf ("=========================================\n");
    return failed == 0 ? 0 : 1;
    }
