/*
 * External allocator-compatibility shim for the curated 1999 KMC model.
 *
 * The archived implementation reads both before and after allocated regions
 * during initialization (confirmed under AddressSanitizer in rxnlist.c:93 and
 * lattice.c:271). Those undefined accesses happened not to abort in the two
 * historical environments but crash reproducibly under modern glibc. This
 * linker wrapper preserves the curated source bytes while supplying explicit,
 * aligned prefix and suffix slack around every malloc allocation. The harness
 * records this shim's hash and the source audit; the wart is not hidden.
 */

#include <stdint.h>
#include <stdlib.h>

#define HISTORICAL_MALLOC_SLACK_BYTES 64U

void *__real_malloc(size_t size);
void __real_free(void *pointer);

void *__wrap_malloc(size_t size)
{
    unsigned char *base;
    if (size > SIZE_MAX - (2U * HISTORICAL_MALLOC_SLACK_BYTES))
        return NULL;
    base = (unsigned char *) __real_malloc(
        size + (2U * HISTORICAL_MALLOC_SLACK_BYTES));
    if (base == NULL)
        return NULL;
    return base + HISTORICAL_MALLOC_SLACK_BYTES;
}

void __wrap_free(void *pointer)
{
    if (pointer != NULL)
        __real_free((unsigned char *) pointer - HISTORICAL_MALLOC_SLACK_BYTES);
}
