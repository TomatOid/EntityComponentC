#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>

unsigned int hammingDist(unsigned int p, unsigned int q)
{
#ifdef __GNUC__
    return __builtin_popcount(p ^ q);
#else
    p ^= q;
    for (q = 0; p; p >>= 1)
        q += p & 1;
    return q;
#endif
}

struct SetNode {
    unsigned int value;
    unsigned int rank;
    struct SetNode *parent;
};

// find the parent-most node connected to x
struct SetNode *find(struct SetNode *x)
{
    struct SetNode *root = x;
    while (root->parent)
        root = root->parent;

    while (x->parent) {
        struct SetNode *parent = x->parent;
        x->parent = root;
        x = parent;
    }

    return root;
}

// returns 1 if any work was done
int unionFind(struct SetNode *x, struct SetNode *y)
{
    // make sure we are root of each set
    x = find(x);
    y = find(y);

    // if they are equal, we know that the two sets
    // are already connected, so we can stop here
    if (x == y)
        return 0;

    // make sure x is the highest rank root
    if (x->rank < y->rank) {
        struct SetNode *temp = x;
        x = y;
        y = temp;
    }

    // make x and y union 
    y->parent = x;
    x->rank += (x->rank == y->rank);
    return 1;
}

struct GraphEdge {
    struct SetNode *verticies[2];
    unsigned char ham_dist;
};

int edgeCompareFn(const void *a, const void *b)
{
    return ((struct GraphEdge *)a)->ham_dist - ((struct GraphEdge *)b)->ham_dist;
}

struct GraphEdge *getMinSpanningForest(unsigned int *numbers, size_t count, struct SetNode **disjoint_set)
{
    // we want to solve this for the complete graph, so we will have
    // count choose 2 edges, so we allocate an array of size count * (count - 1) / 2

    size_t edges_count = count * (count - 1) / 2;
    struct GraphEdge *edges = malloc(sizeof(struct GraphEdge) * edges_count);
    if (!edges)
        return NULL;

    // also allocate enough space for the disjoint-set, and avoid memory leak
    // by setting value of disjoint_set
    if (!disjoint_set)
        return NULL;
    *disjoint_set = malloc(sizeof(struct SetNode) * count);
    struct SetNode *dis_set = *disjoint_set;
    if (!dis_set)
        return NULL;

    // now initialise the disjoint-set
    for (size_t i = 0; i < count; i++)
        dis_set[i] = (struct SetNode) { .rank = 0, .parent = NULL, .value = numbers[i] };

    // fill the edges array with all the pairs of verticies
    size_t k = 0;
    for (size_t i = 0; i < count - 1; i++) {
        for (size_t j = i + 1; j < count; j++, k++) {
            edges[k] = (struct GraphEdge) { .verticies = { dis_set + i, dis_set + j },
                .ham_dist = (unsigned char)hammingDist(numbers[i], numbers[j]) };
        }
    }

    // sort ascending based on length
    qsort(edges, edges_count, sizeof(struct GraphEdge), edgeCompareFn);

    // the minimum spanning tree can have at most count - 1 elements
    struct GraphEdge *result = malloc(sizeof(struct GraphEdge) * (count - 1));
    if (!result)
        return NULL;

    // and now after all that setup, stand by for the algorithm!
    size_t return_count = 0;
    for (size_t i = 0; i < edges_count; i++) { 
        if (unionFind(edges[i].verticies[0], edges[i].verticies[1]))
            result[return_count++] = edges[i];
    }

    free(edges);
    return result;
} 

int main()
{
    unsigned int array[] = { 1, 2, 3, 5, 6, 7 };

    struct SetNode *set = NULL;
    struct GraphEdge *min_tree = getMinSpanningForest(array, sizeof(array) / sizeof(*array), &set);
    
    for (int i = 0; i < sizeof(array) / sizeof(*array) - 1; i++)
        printf("<%d, %d>\n", min_tree[i].verticies[0]->value, min_tree[i].verticies[1]->value);
}
