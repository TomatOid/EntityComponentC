#pragma once
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include "BlockAllocate.h"

typedef struct HashItem
{
    uint64_t key;
    uint64_t value;
    struct HashItem *next;
} HashItem;

// A seprate chaining hash table struct
typedef struct HashTable
{
    HashItem **items;
    HashItem *last;
    BlockPage page;
    size_t len;
    size_t num;
} HashTable;

bool insertToTable(HashTable *table, uint64_t key, uint64_t value)
{
    size_t hash = key % table->len;
    // using a block allocator for this instead
    // of calloc to avoid heap fragmentation
    HashItem *curr = blockAlloc(&table->page);
    if (!curr) return false;
    curr->key = key;
    curr->value = value;
    // insert the item at the beginning of the linked list
    curr->next = table->items[hash];
    table->items[hash] = curr;
    table->num++;
    return true;
}

int findInTable(HashTable *table, uint64_t key, uint64_t *value)
{
    size_t hash = key % table->len;
    HashItem *curr = table->items[hash];

    while (curr && curr->key != key) { curr = curr->next; }

    if (curr) { table->last = curr; *value = curr->value; return 1; }
    else { return 0; }
}

int removeFromTable(HashTable *table, uint64_t key, uint64_t *value)
{
    size_t hash = key % table->len;
    HashItem *curr = table->items[hash];
    HashItem *prev = NULL;
    // loop untill either we hit the end of the list or the keys match
    while (curr && curr->key != key) { prev = curr; curr = curr->next; }

    if (curr) 
    { 
        HashItem *after = curr->next;
        *value = curr->value;
        if (prev) prev->next = after;
        if (!blockFree(&table->page, curr)) return 0; // TODO: use a different value for errors?
        if (!prev) table->items[hash] = after;
        table->num--;
        return 1;
    }
    else return 0;
}

int removeFromTableByValue(HashTable *table, uint64_t key, uint64_t value)
{
    size_t hash = key % table->len;
    HashItem *curr = table->items[hash];
    HashItem *prev = NULL;
    int key_matches_count = 0;
    // loop untill either we hit the end of the list or the keys match
    while (curr && (curr->key != key || curr->value != value)) { key_matches_count += (curr->key == key); prev = curr; curr = curr->next; }
    HashItem *remaining = curr;
    while (remaining) { key_matches_count += remaining->key == key; remaining = remaining->next; }

    if (curr) 
    { 
        HashItem *after = curr->next;
        if (prev) prev->next = after;
        if (!blockFree(&table->page, curr)) return 0;
        if (!prev) table->items[hash] = after;
        table->num--;
        return key_matches_count;
    } else return 0;
}
