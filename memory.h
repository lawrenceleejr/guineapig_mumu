struct struct_memory_pointer
{
  void *memory;
  struct struct_memory_pointer *next;
};

typedef struct struct_memory_pointer MEMORY_POINTER;

struct struct_memory_account
{
  long amount,calls,size;
  MEMORY_POINTER *first;
};

typedef struct struct_memory_account MEMORY_ACCOUNT;

extern int store_memory_monitor[10];

void init_memory_account(MEMORY_ACCOUNT*,long);
void* get_memory(MEMORY_ACCOUNT*,int);
void clear_memory_account(MEMORY_ACCOUNT*);
long memory_account_amount(MEMORY_ACCOUNT*);
long memory_account_size(MEMORY_ACCOUNT*);
long memory_account_calls(MEMORY_ACCOUNT*);
extern MEMORY_ACCOUNT m_account;
