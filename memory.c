#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**********************************************************************/
/* memory_account routines                                            */
/**********************************************************************/

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

/* Storage and routines for the memory requirements */

int store_memory_monitor[10];

void init_memory_account(MEMORY_ACCOUNT *account,long size)
{
  account->size=size;
  account->calls=0;
  account->amount=0;
  account->first=NULL;
}

void* get_memory(MEMORY_ACCOUNT *account,int size)
{
  MEMORY_POINTER *temp;
  if (account->size>0)
    {
      if (account->size<account->amount+size)
	{
	  fprintf(stderr,"error: accout exceeded\n");
	  fflush(stderr);
	  exit(10);
	}
    }
  account->calls++;
  account->amount += size;
  temp=(MEMORY_POINTER*)malloc(sizeof(MEMORY_POINTER));
  if (temp==NULL){
      fprintf(stderr,"not eneough memory\n");
      fprintf(stderr,"%d bytes required\n",size);
      fflush(stderr);
      exit(10);
  }
  temp->next=account->first;
  temp->memory=malloc(size);
  if(temp->memory==NULL){
      fprintf(stderr,"not eneough memory\n");
      fprintf(stderr,"%d bytes required\n",size);
      fflush(stderr);
      exit(10);
  }
  account->first=temp;
  return temp->memory;
}

void clear_memory_account(MEMORY_ACCOUNT *account)
{
  MEMORY_POINTER *temp;
  account->amount=0;
  account->calls=0;
  temp=account->first;
  while(temp!=NULL)
    {
      account->first=temp->next;
      free(temp->memory);
      free(temp);
      temp=account->first;
    }
}

long memory_account_amount(MEMORY_ACCOUNT *account)
{
  return account->amount;
}

long memory_account_size(MEMORY_ACCOUNT *account)
{
  return account->size;
}

long memory_account_calls(MEMORY_ACCOUNT *account)
{
  return account->calls;
}

MEMORY_ACCOUNT m_account;
