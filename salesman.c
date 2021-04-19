#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include "HashTable.h"

uint64_t hammingDist(uint64_t p, uint64_t q)
{
#ifdef __GNUC__
    return __builtin_popcountl(p ^ q);
#else
    p ^= q;
    for (q = 0; p; p >>= 1)
        q += p & 1;
    return q;
#endif
}

struct SetNode {
    uint64_t value;
    size_t rank;
    unsigned char degree;
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
int unionFind(struct SetNode *x_base, struct SetNode *y_base)
{
    // make sure we are root of each set
    struct SetNode *x = find(x_base);
    struct SetNode *y = find(y_base);

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
    x_base->degree++;
    y_base->degree++;
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

struct GraphEdge *getMinSpanningForest(uint64_t *numbers, size_t count, struct SetNode **disjoint_set)
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

    // sort ascending based on hamming distance
    qsort(edges, edges_count, sizeof(struct GraphEdge), edgeCompareFn);

    // the minimum spanning tree can have at most count - 1 elements
    struct GraphEdge *result = malloc(sizeof(struct GraphEdge) * (count - 1));
    if (!result)
        return NULL;

    // and now after all that setup, stand by for the algorithm!
    size_t return_count = 0;
    for (size_t i = 0; i < edges_count; i++) { 
        if (unionFind(edges[i].verticies[0], edges[i].verticies[1])) {
            result[return_count++] = edges[i];
        }
    }

    free(edges);
    return result;
} 

// Some reasoning behind this algorithm:
//
// It seems like most perfect matching algorithms think about the graph
// additively, meaning that they start with an empty graph and add connections
// in such a way that will: 
//  (a) increase the total weight as little as possible
//  (b) make the degree of each vertex 1
//
// It would be easier to think of this algorithm as subtractive, meaning that
// you start with a complete graph and eliminate edges in such a way that:
//  (a) isolates a set of verticies from the rest of the graph (making them
//  degree 1)
//  (b) decreases the total weight of the graph as much as possible
//
// This step can be applied recursively to the larger subgraph to produce the same result
// as an additive algorithm
int getMinMatching(uint64_t *odd_nodes, size_t count)
{
    if (!count || count & 1) 
        return 0;

    uint64_t distance_sums[count];

    // now precompute the distances, this reduces the time complexity from O(n^4) to O(n^3)
    for (size_t i = 0; i < count; i++)
        for (size_t j = distance_sums[i] = 0; j < count; j++)
            distance_sums[i] += hammingDist(odd_nodes[i], odd_nodes[j]);

    for (; count >= 2; count -= 2) {
        ssize_t min_oppertunity_cost = SSIZE_MAX;
        size_t min_i, min_j = 0;
        // loop over all edges
        for (size_t i = 0; i < count - 1; i++) {
            for (size_t j = i + 1; j < count; j++) {
                // calculate the oppertunity cost of connecting these two nodes
                ssize_t cost = 2 * hammingDist(odd_nodes[i], odd_nodes[j]) - distance_sums[i] - distance_sums[j];
                
                if (cost < min_oppertunity_cost) {
                    min_oppertunity_cost = cost;
                    min_i = i;
                    min_j = j;
                }
            }
        }
        // update distance_sums to account for the removal
        for (size_t i = 0; i < count; i++)
            distance_sums[i] -= hammingDist(odd_nodes[i], odd_nodes[min_i]) + hammingDist(odd_nodes[i], odd_nodes[min_j]);
        distance_sums[min_j] = distance_sums[count - 1];
        distance_sums[min_i] = distance_sums[count - 2];
        // swap the min i and j to the end
        uint64_t t = odd_nodes[min_j];
        odd_nodes[min_j] = odd_nodes[count - 1];
        odd_nodes[count - 1] = t;
        t = odd_nodes[min_i];
        odd_nodes[min_i] = odd_nodes[count - 2];
        odd_nodes[count - 2] = t;
    }
    return 1;
}

uint64_t *getOddDegree(struct SetNode *set, size_t count, size_t *odd_degree_count)
{
    // first count the number of odd degree elements
    size_t odd_count = 0;
    for (size_t i = 0; i < count; i++)
        odd_count += set[i].degree & 1;

    *odd_degree_count = odd_count;

    uint64_t *result = malloc(sizeof(uint64_t) * odd_count);
    if (!result)
        return NULL;
    
    odd_count = 0;
    for (size_t i = 0; i < count; i++)
        if (set[i].degree & 1)
            result[odd_count++] = set[i].value;

    return result;
}

uint64_t *getEulerPath(uint64_t *edges, size_t count)
{
    uint64_t *solution = NULL;
    // first set up an ajacency hash table
    BlockPage page;
    if (!makePage(&page, count, sizeof(HashItem)))
        return NULL;
    HashTable adj_table = (HashTable) { .items = calloc(count, sizeof(HashItem *)),
        .page = page, .len = count };

    if (!adj_table.items)
        goto error0;

    for (size_t i = 0; i < count; i += 2) {
        insertToTable(&adj_table, edges[i], edges + i + 1);
        insertToTable(&adj_table, edges[i + 1], edges + i);
    }

    // allocate space for the stacks
    ssize_t stack_top = 0;
    uint64_t **stack = malloc(sizeof(uint64_t *) * count);
    if (!stack)
        goto error1;
    size_t solution_top = 0;
    solution = malloc(sizeof(uint64_t) * count);
    if (!solution)
        goto error2;
    
    uint64_t *current_node = edges;
    do {
        uint64_t *next_node = removeFromTable(&adj_table, *current_node);
        if (next_node) {
            removeFromTableByValue(&adj_table, *next_node, current_node);
            stack[stack_top++] = current_node;
            current_node = next_node;
        }
        else { 
            solution[solution_top++] = *current_node;
            current_node = stack[--stack_top];
        }
    } while (stack_top > 0);

    printf("solution_top = %lu, count = %lu\n", solution_top, count);

    // think i know why people like defer now
    // TODO: maybe do one giant alloc instead of a bunch of little ones,
    // although this will probably increase the chance of failure in a
    // fragmented heap
error2:
    free(stack);
error1:
    free(adj_table.items);
error0:
    free(page.pool);
    free(page.free);
    return solution;
}

#define array_size(arr) (sizeof(arr) / sizeof(*arr))

int main()
{
    //uint64_t array[] = { 1, 2, 3, 5, 6, 7, 12, 13, 16, 17, 19, 22, 31 };
    //uint64_t array[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    uint64_t array[8];
    for (uint64_t i = 0; i < array_size(array); i += 1)
        array[i] = i;

    struct SetNode *set = NULL;
    struct GraphEdge *min_tree = getMinSpanningForest(array, array_size(array), &set);

    for (int i = 0; i < array_size(array) - 1; i++)
        printf("<%lu, %lu>\n", min_tree[i].verticies[0]->value, min_tree[i].verticies[1]->value);

    //for (int i = 0; i < array_size(array); i++)
    //    printf("value: %lu, degree: %d\n", set[i].value, set[i].degree);

    size_t odd_count = 0;
    uint64_t *odds = getOddDegree(set, sizeof(array) / sizeof(*array), &odd_count);
    
    getMinMatching(odds, odd_count);

    for (int i = 0; i + 1 < odd_count; i += 2)
        printf("[%lu, %lu]\n", odds[i], odds[i + 1]);

    // append the two arrays together
    // TODO: correct spaget
    size_t total_len = 2 * (array_size(array) - 1) + odd_count;
    uint64_t *euler_in = malloc(sizeof(uint64_t) * total_len);
    int index = 0;
    for (; index / 2 < array_size(array) - 1; index += 2) {
        euler_in[index] = min_tree[index / 2].verticies[0]->value;
        euler_in[index + 1] = min_tree[index / 2].verticies[1]->value;
    }
    for (int i = 0; index < total_len; i++, index++) {
        euler_in[index] = odds[i];
    }
    
    uint64_t *euler = getEulerPath(euler_in, total_len);
    for (int i = 0; i < total_len - 1; i++) {
        printf("%lu, ", euler[i]);
    }

}
