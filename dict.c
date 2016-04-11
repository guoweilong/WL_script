// dict.h
//
// Copyright (c) 2015 c4605 <bolasblack@gmail.com>
// Code from : https://github.com/bolasblack/dict.c/tree/master/src


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define DICT_MALLOC malloc
#define DICT_FREE free
// dict_t pair struct.
typedef struct dict_pair {
  struct dict_pair *prev;
  struct dict_pair *next;
  char *key;
  void *val;
} dict_pair_t;
// dict_t struct.
typedef struct {
  dict_pair_t *head;
  dict_pair_t *tail;
  void (*free)(char *key, void *val);
} dict_t;
// dict_t iterator struct.
typedef struct {
  dict_pair_t *next;
} dict_iterator_t;
// Node prototypes.
dict_pair_t *
dict_pair_new(char *key, void *val);
// dict_t prototypes.
dict_t * dict_new();
dict_pair_t * dict_set(dict_t *self, char *key, void *val);
dict_pair_t * dict_get(dict_t *self, char *key);
void dict_remove(dict_t *self, char *key);
void dict_destroy(dict_t *self);
// dict_t iterator prototypes.
dict_iterator_t * dict_iterator_new(dict_t *dict);
dict_pair_t * dict_iterator_next();
void dict_iterator_destroy(dict_iterator_t *self);
// Allocate a new dict_iterator_t. NULL on failure.
dict_iterator_t * dict_iterator_new(dict_t *dict) {
  dict_pair_t *pair = dict->head;
  dict_iterator_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_iterator_t))))
    return NULL;
  self->next = pair;
  return self;
}
// Return the next dict_pair_t or NULL when no more
// nodes remain in the dict.
dict_pair_t * dict_iterator_next(dict_iterator_t *self) {
  dict_pair_t *curr = self->next;
  if (curr)
    self->next = curr->next;
  return curr;
}
// Free the dict iterator.
void dict_iterator_destroy(dict_iterator_t *self) {
  DICT_FREE(self);
  self = NULL;
}
// dict.c
// Allocates a new dict_pair_t. NULL on failure.
dict_pair_t * dict_pair_new(char *key, void *val) {
  dict_pair_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_pair_t))))
    return NULL;
  self->prev = NULL;
  self->next = NULL;
  self->key = key;
  self->val = val;
  return self;
}
// Allocate a new dict_t. NULL on failure.
dict_t * dict_new() {
  dict_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_t))))
    return NULL;
  self->head = NULL;
  self->tail = NULL;
  self->free = NULL;
  return self;
}
// Free the dict.
void dict_destroy(dict_t *self) {
  dict_pair_t *next;
  dict_pair_t *curr = self->head;
  while (curr) {
    next = curr->next;
    if (self->free) self->free(curr->key, curr->val);
    DICT_FREE(curr);
    curr = next;
  }
  DICT_FREE(self);
}
// Add a key-value pair to dict.
dict_pair_t * dict_set(dict_t *self, char *key, void *val) {
  dict_pair_t *pair = dict_get(self, key);
  if (pair) {
    if (self->free) self->free(pair->key, pair->val);
    pair->val = val;
  } else {
    pair = dict_pair_new(key, val);
    if (self->head) {
      pair->prev = self->tail;
      pair->next = NULL;
      self->tail->next = pair;
      self->tail = pair;
    } else {
      self->head = self->tail = pair;
      pair->prev = pair->next = NULL;
    }
  }
  return pair;
}
// Return the node associated to val or NULL.
dict_pair_t * dict_get(dict_t *self, char *key) {
  dict_iterator_t *it = dict_iterator_new(self);
  dict_pair_t *pair;
  while ((pair = dict_iterator_next(it))) {
    if (strcmp(key, pair->key) == 0) {
      dict_iterator_destroy(it);
      return pair;
    }
  }
  dict_iterator_destroy(it);
  return NULL;
}
// Remove the given node from the dict, freeing it and it's value.
void dict_remove(dict_t *self, char *key) {
  dict_pair_t *pair = dict_get(self, key);
  if (!pair) return;
  pair->prev
    ? (pair->prev->next = pair->next)
    : (self->head = pair->next);
  pair->next
    ? (pair->next->prev = pair->prev)
    : (self->tail = pair->prev);
  if (self->free) self->free(pair->key, pair->val);
  DICT_FREE(pair);
}



int main(){
    dict_t * Dict = dict_new();
    dict_set(Dict, "A", "1");
    dict_set(Dict, "B", "1");
    dict_set(Dict, "C", "1");
    char key[100] = "A";
    dict_pair_t * DictPair = dict_get(Dict, key);
    printf("Dict[%s]=%s\n", DictPair->key, (char *) DictPair->val);
	dict_destroy(Dict);
    dict_set(Dict, "A", "2");
    printf("Dict[%s]=%s\n", DictPair->key, (char *) DictPair->val);
    return 1;
}


