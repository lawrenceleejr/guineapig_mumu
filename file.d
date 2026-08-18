#include <ctype.h>

/* datei handling */

#define max_buffer 512
#define max_buffer2 10000

/*
#define GETCH getc(stdin)
#define UNGETCH(a) ungetc(a,stdin)
*/
#define GETCH get_char()
#define UNGETCH(a) unget_char(a)

#define CONTENTS(a) (((a).type==T_DOUBLE)?(a).value.d:(a).value.i)

struct {char *buffer;} input;

int get_char()
{
  int tmp;
  tmp=*(input.buffer);
  input.buffer++;
  return (int)tmp;
}

void
unget_char(char a)
{
  input.buffer--;
}

void
store_buffer(char *buffer)
{
  input.buffer=buffer;
}

typedef enum {T_DOUBLE,T_INT,T_STRING,T_DOUBLE_ARRAY,T_INT_ARRAY,
	      T_DOUBLE_MIRROR,T_INT_MIRROR,T_DOUBLE_2,T_INT_2} VAR_TYPE;

typedef struct {VAR_TYPE type;union {double d;long int i;char *s;} value;}
VALUE;

int
value_to_int(VALUE *val)
{
	if(val->type!=T_INT){
	  fprintf(stderr,"value is not an integer\n");
	  exit(1);
	}
	return (val->value).i;
}

double
value_to_double(VALUE *val)
{
	if(val->type!=T_DOUBLE){
	  fprintf(stderr,"value is not a double\n");
	  exit(1);
	}
	return (val->value).d;
}

void
int_to_value(int i,VALUE *val)
{
	val->type=T_INT;
	(val->value).i=i;
}

void
double_to_value(double d,VALUE *val)
{
	val->type=T_DOUBLE;
	(val->value).d=d;
}

typedef struct {
char name[32];
VAR_TYPE type;
union {double d;
       long int i;
       double *dp;
       long int *ip;
       double dm[2];
       long int im[2];
       char *string;} value;
} VARIABLE;

typedef int INDEX;

typedef enum {START,MULT,DIV,PLUS,MINUS,END,END_INPUT,ERROR,ASSIGN,LB,RB,VAL,
	      VAR_NAME,STRING} TOKEN;

struct{int number,max_number; VARIABLE *heap;} var_heap;

void read_error(char error[])
{
    fprintf(stderr,"Error:\n%s\n",error);
}

void
read_value_unit(char *buffer)
{
  int letter;
  letter=GETCH;
  while(isdigit(letter)){
    *(buffer++)=letter;
    letter=GETCH;
  }
  if (letter=='.'){
    *(buffer++)=letter;
    letter=GETCH;
    while(isdigit(letter)){
      *(buffer++)=letter;
      letter=GETCH;
    }
  }
  if (letter=='e'||letter=='E'){
    *(buffer++)=letter;
    letter=GETCH;
    if (letter=='-'||letter=='+'){
      *(buffer++)=letter;
      letter=GETCH;
    }
    while(isdigit(letter)){
      *(buffer++)=letter;
      letter=GETCH;
    }
  }
  *buffer=(char)0;
  UNGETCH(letter);
}

int convert_to_value(char buffer[],VALUE *val)
{
    int letter;
    char *end;
    int i=0,n=0,isdouble=0;

    letter=buffer[0];
    if (!isdigit(letter)&&letter!='.'&&letter!='-'){
	read_error("Value is not a number");
	return 1;
    }
    while(buffer[++i]!='\0'){
	if (buffer[i]=='.'||buffer[i]=='e'||buffer[i]=='E'){
	    isdouble=1;
	}
    }
    if (isdouble){
	val->type=T_DOUBLE;
	val->value.d=strtod(buffer,&end);
    }
    else{
	val->type=T_INT;
	val->value.i=strtol(buffer,&end,10);
    }
    if (*end!='\0') return 1;
    return 0;
}

TOKEN read_value(VALUE *val)
{
    int i=0,n=1000;
    char buffer[1000];
    int letter;

    read_value_unit(buffer);
    if (convert_to_value(buffer,val)) return ERROR;
    return VAL;
/*
    letter=GETCH;
    while(letter!=EOF){
	if (!isdigit(letter)&&letter!='.'&&letter!='e'&&letter!='E'){
	    UNGETCH(letter);
	    buffer[i]='\0';
	    if (convert_to_value(buffer,val)) return ERROR;
	    return VAL;
	}
	buffer[i++]=letter;
	letter=GETCH;   
	if (i==n) {
	    return ERROR;
	}
    }
    return END_INPUT;
*/
}

TOKEN read_string(char *buffer,int n)
{
    int letter;
    int i;

    letter=GETCH;
    if (letter!=(int)'"') {
	read_error("String definition starts with wrong letter");
	return ERROR;
    }
    for (i=1;i<n;i++){
	letter=GETCH;
	if (letter=='"'){
	    return STRING;
	}
	if (letter==EOF){
	    read_error("File ended while scanning of string");
	    return ERROR;
	}
	*buffer++=letter;
    }
    read_error("String too long");
    return ERROR;
}

void
def_variable(VARIABLE *var,VAR_TYPE type)
{
    var->type=type;
}

void
set_variable(VARIABLE *var,VALUE val)
{
    switch(var->type){
    case T_DOUBLE:
	var->value.d=CONTENTS(val);
	break;
    case T_INT:
	var->value.i=CONTENTS(val);
	break;
    case T_DOUBLE_MIRROR:
	var->value.dm[0]=CONTENTS(val);
	var->value.dm[1]=-CONTENTS(val);
	break;
    case T_INT_MIRROR:
	var->value.im[0]=CONTENTS(val);
	var->value.im[1]=-CONTENTS(val);
	break;
    case T_DOUBLE_2:
	var->value.dm[0]=CONTENTS(val);
	var->value.dm[1]=CONTENTS(val);
	break;
    case T_INT_2:
	var->value.im[0]=CONTENTS(val);
	var->value.im[1]=CONTENTS(val);
	break;
    }
}

void
set_variable_string(VARIABLE *var,char cont[])
{
    int n=0,i;
    if (var->type!=T_STRING){
	fprintf(stderr,"Variable is not of type string\n");
	exit(1);
    }

    while(cont[n++]!='\0') ;

    var->value.string=(char*)get_memory(&m_account,sizeof(char)*n);

    strcpy(var->value.string,cont);
}

void
get_variable(VARIABLE *var,VALUE *val)
{
    switch(var->type){
    case T_DOUBLE:
	val->value.d=var->value.d;
	val->type=var->type;
	break;
    case T_INT:
	val->value.i=var->value.i;
	val->type=var->type;
	break;
    case T_DOUBLE_MIRROR:
    case T_INT_MIRROR:
    case T_DOUBLE_2:
    case T_INT_2:
	fprintf(stderr,"Error: Variable <%s> not a single value\n",var);
	fprintf(stderr,"Error occured in get_variable\n");
	exit(1);
    case T_STRING:
	val->type=T_STRING;
	val->value.s=var->value.string;
    }
}

void
get_variable_element(VARIABLE *var,INDEX index,VALUE *val)
{
    switch(var->type){
    case T_DOUBLE:
    case T_INT:
	fprintf(stderr,"Error: Variable <%s> is a single value\n",var);
	fprintf(stderr,"Error occured in get_variable_element\n");
	exit(1);
    case T_DOUBLE_MIRROR:
    case T_DOUBLE_2:
	val->type=T_DOUBLE;
	if (index==1){
	    val->value.d=var->value.dm[0];
	}
	else{
	    if (index==2){
		val->value.d=var->value.dm[1];
	    }
	    else{
	        fprintf(stderr,"Error: Variable <%s> has only elements 1 \
and 2\n",var);
		fprintf(stderr,"Error occured in get_variable_element\n");
		exit(1);
	    }
	}
	break;
    case T_INT_MIRROR:
    case T_INT_2:
	val->type=T_INT;
	if (index==1){
	    val->value.i=var->value.im[0];
	}
	else{
	    if (index==2){
		val->value.i=var->value.im[1];
	    }
	    else{
	        fprintf(stderr,"Error: Variable <%s> has only elements 1 \
and 2\n",var);
		fprintf(stderr,"Error occured in get_variable_element\n");
		exit(1);
	    }
	}
    }
}

void
set_variable_element(VARIABLE *var,INDEX index,VALUE val)
{
    char msg[200];

    switch(var->type){
    case T_DOUBLE:
    case T_INT:
        fprintf(stderr,"Error: Variable <%s> is single value\n",var);
	fprintf(stderr,"Error occured in set_variable_element\n");
	exit(1);
	break;
    case T_DOUBLE_MIRROR:
    case T_DOUBLE_2:
	if (index==1){
	    var->value.dm[0]=CONTENTS(val);
	}
	else{
	    if (index==2){
		var->value.dm[1]=CONTENTS(val);
	    }
	    else{
                sprintf(msg,"Error: Variable <%s> contains only elements 1 \
and 2\n",var);
                read_error(msg);
        	fprintf(stderr,"Error occured in set_variable_element\n");
		exit(1);
	    }
	}
	break;
    case T_INT_MIRROR:
    case T_INT_2:
	if (index==1){
	    var->value.im[0]=CONTENTS(val);
	}
	else{
	    if(index==2){
		var->value.im[1]=CONTENTS(val);
	    }
	    else{
                sprintf(msg,"Error: Variable <%s> contains only elements 1 \
and 2\n",var);
                read_error(msg);
        	fprintf(stderr,"Error occured in set_variable_element\n");
		exit(1);
	    }
	}
	break;
    }
}

/*
TOKEN read_skip()
{
    int letter;

    while(isspace(letter=GETCH)) ;
    if (letter==EOF) return END_INPUT;
    UNGETCH(letter);
}
*/

void read_skip()
{
    int letter;

    while(isspace(letter=GETCH)) ;
    if (letter==EOF) return;
    UNGETCH(letter);
}


void
print_variable(VARIABLE *var)
{
    switch (var->type){
    case T_DOUBLE:
	printf("DOUBLE:%g\n",var->value.d);
	break;
    case T_INT:
	printf("INT: %d\n",var->value.i);
	break;
    case T_DOUBLE_2:
    case T_DOUBLE_MIRROR:
	printf("DOUBLE.1: %g\n",var->value.dm[0]);
	printf("DOUBLE.2: %g\n",var->value.dm[1]);
	break;
    case T_INT_2:
    case T_INT_MIRROR:
	printf("INT.1: %d\n",var->value.im[0]);
	printf("INT.2: %d\n",var->value.im[1]);
	break;
    case T_STRING:
	printf("STRING: %s\n",var->value.string);
	break;
    }
}

VARIABLE* find_named_variable(char name[])
{
    int i;
    for (i=0;i<var_heap.number;i++){
	if (!strcmp(name,var_heap.heap[i].name)){
	    return &(var_heap.heap[i]);
	}
    }
    return NULL;
}

void
def_named_variable(char name[],VAR_TYPE type)
{
    if (find_named_variable(name)==NULL){
	if (var_heap.max_number-var_heap.number>0){
	    strncpy((var_heap.heap[var_heap.number]).name,name,32);
	    (var_heap.heap[var_heap.number]).type=type;
	    var_heap.number++;
	}
	else{
	    fprintf(stderr,"Too many variables defined\n");
	    exit(1);
	}
    }
    else{
	fprintf(stderr,"Variable <%s> already exists\n",name);
	exit(1);
    }
}

void
print_named_variable(char name[])
{
    VARIABLE *var;
    var=find_named_variable(name);
    if (var!=NULL){
	printf("%s\n",name);
	print_variable(var);
    }
    else{
	fprintf(stderr,"Variable <%s> does not exist\n",name);	
    }
}

void
init_named_variable(int max_entry)
{
    var_heap.number=0;
    var_heap.max_number=max_entry;
    var_heap.heap=(VARIABLE*)get_memory(&m_account,sizeof(VARIABLE)*max_entry);
}

void
set_named_variable(char name[],VALUE val)
{
    char buffer[1000],*point,*end;
    int is_simple;
    INDEX index;
    VARIABLE *var;

    strncpy(buffer,name,1000);
    point=buffer;
    while(isalnum(*point)||*point=='_') point++;
    if (*point=='\0') {
	is_simple=1;
    }
    else{
	if(*point=='.'){
	    is_simple=0;
	}
	else{
	    fprintf(stderr,"Error: Variable <%s> is only single value\n",name);
            fprintf(stderr,"Error occured in set_named_variable\n");
	    exit(1);
	}
    }
    *point='\0';
    var=find_named_variable(buffer);
    if(var==NULL) {
	fprintf(stderr,"Variable <%s> is not defined\n",buffer);
	exit(1);
    }
    if (is_simple){
	set_variable(var,val);
    }
    else{
	point++;
	index=strtol(point,&end,10);
	set_variable_element(var,index,val);
    }
}

void
set_named_variable_string(char name[],char cont[])
{
    char buffer[1000],*point,*end;
    int is_simple;
    INDEX index;
    VARIABLE *var;

    strncpy(buffer,name,1000);
    point=buffer;
    while(isalnum(*point)||*point=='_') point++;
    if (*point=='\0') {
	is_simple=1;
    }
    else{
	fprintf(stderr,"Variable <%s> is not of type string",name);
	exit(1);
    }
    *point='\0';
    var=find_named_variable(buffer);
    if(var==NULL) {
	fprintf(stderr,"Variable <%s> is not defined\n",buffer);
	exit(1);
    }
    if (var->type!=T_STRING){
	fprintf(stderr,"Variable <%s> is not of type string",name);
	exit(1);
    }
    set_variable_string(var,cont);
}

int get_named_variable(char name[],VALUE *val)
{
    char buffer[1000],*point,*end;
    int is_simple;
    INDEX index;
    VARIABLE *var;

    strncpy(buffer,name,1000);
    point=buffer;
    while(isalnum(*point)||*point=='_') point++;
    if (*point=='\0') {
	is_simple=1;
    }
    else{
	if(*point=='.'){
	    is_simple=0;
	}
	else{
	    is_simple=1;
	}
    }
    *point='\0';
    var=find_named_variable(buffer);
    if (var==NULL) {
/* scd temp */
        printf("Variable <%s> not found\n",buffer);
        exit(1);
	return 1;
    }
    if (is_simple){
	get_variable(var,val);
    }
    else{
	point++;
	index=strtol(point,&end,10);
	get_variable_element(var,index,val);
    }
    return 0;
}

TOKEN read_name0(char *buffer,int n)
{
    int i=0,letter,isvar_element=0;

    n--;
    letter=GETCH;
    if (!isalpha(letter)){
	return ERROR;
    }
    buffer[i++]=letter;
    while((letter=GETCH)!=EOF){
	if (i==n){
	    fprintf(stderr,"Name too long in read_name\n");
	    exit(1);
	}
	if (isalnum(letter)||letter=='_'||letter=='.'){
	    buffer[i++]=letter;
	}
	else{
	    UNGETCH(letter);
	    buffer[i]='\0';
	    return VAR_NAME;
	}
    }
    return END_INPUT;
}

TOKEN read_name(char *buffer,int n)
{
    TOKEN token;
    token=read_name0(buffer,n);
    return token;
}

TOKEN read_token0(VALUE *val)
{
    int letter,i=0,n=1000,tmp;
    char buffer[1000];
    TOKEN token;

    read_skip();
    letter=GETCH;
    buffer[0]='\0';
    if(letter==EOF){
	return END_INPUT;
    }
    if (isalpha(letter)){
	UNGETCH(letter);
	token=read_name(buffer,n);
	tmp=get_named_variable(buffer,val);
	if(tmp) return ERROR;
	return VAL;
    }
    if (isdigit(letter)||letter=='.'){
	UNGETCH(letter);
	return read_value(val);
    }
    switch(letter){
    case '{':
	return START;
    case '*':
	return MULT;
    case '/':
	return DIV;
    case '+':
	return PLUS;
    case '-':
	return MINUS;
    case '}':
	return END_INPUT;
    case ';':
	return END;
    case '=':
	return ASSIGN;
    case '"':
	UNGETCH(letter);
	return read_string(buffer,n);
    }
    return ERROR;
}

void
print_token(TOKEN token)
{
    printf("<");
    switch(token){
    case ERROR:
	printf("ERROR");
	break;
    case STRING:
	printf("STRING");
	break;
    case RB:
	printf("RB");
	break;
    case LB:
	printf("LB");
	break;
    case VAL:
	printf("VAL");
	break;
    case VAR_NAME:
	printf("VAR_NAME");
	break;
    case MULT:
	printf("MULT");
	break;
    case DIV:
	printf("DIV");
	break;
    case PLUS:
	printf("PLUS");
	break;
    case MINUS:
	printf("MINUS");
	break;
    case END:
	printf("END");
	break;
    case END_INPUT:
	printf("END_INPUT");
	break;
    case ASSIGN:
	printf("ASSIGN");
	break;
    }
    printf(">");
}

TOKEN read_token(VALUE *val)
{
    TOKEN token;
    token=read_token0(val);
/*    print_token(token);*/
    return token;
}

TOKEN read_token_name0(char buffer[],int n)
{
    int letter,i=0;
    TOKEN token;

    read_skip();
    letter=GETCH;
    buffer[0]='\0';
    if(letter==EOF){
	return END_INPUT;
    }
    if (isalpha(letter)){
	UNGETCH(letter);
	token=read_name(buffer,n);
	return token;
    }
    if (isdigit(letter)||letter=='.'){
	UNGETCH(letter);
	return ERROR;
    }
    switch(letter){
    case '*':
	return MULT;
    case '/':
	return DIV;
    case '+':
	return PLUS;
    case '-':
	return MINUS;
    case '}':
	return END_INPUT;
    case ';':
	return END;
    case '=':
	return ASSIGN;
    case '"':
	UNGETCH(letter);
	return read_string(buffer,n);
    }
    return ERROR;
}

TOKEN read_token_name(char buffer[],int n)
{
    TOKEN token;
    token=read_token_name0(buffer,n);
/*    printf("(1");print_token(token);printf(")2");*/
    return token;
}

TOKEN read_lvalue(char name[],int n)
{
    TOKEN token;
/*printf("[1");*/
    token=read_token_name(name,n);
    if (token==END_INPUT) return END_INPUT;
/*printf("]2");*/
    if (token!=VAR_NAME) return ERROR;
/*    printf ("<%s>\n",name);*/
    return VAR_NAME;
}

void
invert_value(VALUE *val)
{
    switch(val->type){
    case T_DOUBLE:
	val->value.d=-val->value.d;
	break;
    case T_INT:
	val->value.i=-val->value.i;
	break;
    default:
	break;
    }
}

void
add_value(VALUE *val1,VALUE *val2,VALUE *val3)
{
    if (val1->type==T_DOUBLE||val2->type==T_DOUBLE){
	val3->value.d=CONTENTS(*val1)+CONTENTS(*val2);
	val3->type=T_DOUBLE;
    }
    else{
	val3->value.i=CONTENTS(*val1)+CONTENTS(*val2);
	val3->type=T_INT;
    }
}

void
sub_value(VALUE *val1,VALUE *val2,VALUE *val3)
{
    if (val1->type==T_DOUBLE||val2->type==T_DOUBLE){
	val3->value.d=CONTENTS(*val1)-CONTENTS(*val2);
	val3->type=T_DOUBLE;
    }
    else{
	val3->value.i=CONTENTS(*val1)-CONTENTS(*val2);
	val3->type=T_INT;
    }
}

void
mult_value(VALUE *val1,VALUE *val2,VALUE *val3)
{
    if (val1->type==T_DOUBLE||val2->type==T_DOUBLE){
	val3->value.d=CONTENTS(*val1)*CONTENTS(*val2);
	val3->type=T_DOUBLE;
    }
    else{
	val3->value.i=CONTENTS(*val1)*CONTENTS(*val2);
	val3->type=T_INT;
    }
}

void
div_value(VALUE *val1,VALUE *val2,VALUE *val3)
{
    if (val1->type==T_DOUBLE||val2->type==T_DOUBLE){
	val3->value.d=CONTENTS(*val1)/CONTENTS(*val2);
	val3->type=T_DOUBLE;
    }
    else{
	val3->value.i=CONTENTS(*val1)/CONTENTS(*val2);
	val3->type=T_INT;
    }
}

TOKEN read_expression(VALUE*);

TOKEN read_primitive(VALUE *val)
{
    TOKEN token;
    VALUE val1;
/*printf("!1");*/
    token=read_token(val);
    switch(token){
    case VAL:
	break;
    case MINUS:
	token=read_primitive(val);
	invert_value(val);
/*printf("!3");*/
	return token;
    case LB:
	token=read_expression(val);
	if (token!=VAL) return ERROR;
	if ((token=read_token(&val1))!=RB) return ERROR;
	return token;
    default:
	return ERROR;
    }
/*printf("!2");*/
    return read_token(&val1);
}

TOKEN read_summand(VALUE *val)
{
    TOKEN token;
    VALUE val1;
/*printf("(");*/
    token=read_primitive(val);
    while((token==MULT)||(token==DIV)){
	switch(token){
	case MULT:
	    token=read_primitive(&val1);
	    mult_value(val,&val1,val);
	    break;
	case DIV:
	    token=read_primitive(&val1);
	    div_value(val,&val1,val);
	    break;
	default:
	    break;
	}
    }
/*printf(")");*/
    return token;
}

TOKEN read_expression(VALUE *val)
{
    char buffer[1000];
    VALUE val1;
    TOKEN token;
/*printf("[");*/
    token=read_summand(val);
    if ((token==MINUS)||(token==PLUS)){
	if (token==MINUS){
	    token=read_expression(&val1);
	    if (token==ERROR) return ERROR;
	    sub_value(val,&val1,val);
	}
	else{
	    token=read_expression(&val1);
	    if (token==ERROR) return ERROR;
	    add_value(val,&val1,val);
	}
    }
/*printf("]");*/
    return token;
}

TOKEN read_command()
{
    char lval[1000],rval[1000],dummy[2];
    VALUE val;
    VARIABLE *tmp;
    TOKEN token;

    if ((token=read_lvalue(lval,1000))==END_INPUT) return END_INPUT;
    if (token!=VAR_NAME) return ERROR;
    if (read_token(&val)!=ASSIGN) return ERROR;
    tmp=find_named_variable(lval);
    if (tmp!=NULL){
	if (tmp->type==T_STRING){
	    read_string(rval,1000);
	    set_named_variable_string(lval,rval);
	    if (read_token(&val)!=END) return ERROR;
	    return ASSIGN;
	}
    }
    if (read_expression(&val)!=END){
	return ERROR;
    }
    set_named_variable(lval,val);
    return ASSIGN;
}

void
define_values(char *line)
{
  TOKEN token;
  VALUE value;
  store_buffer(line);
  if ((token=read_token(&value))!=START){
	printf("error: no proper start\n");
	exit(1);
  }
  while((token=read_command())!=END_INPUT) {
/*	printf("{");
	print_token(token);
	printf("}");*/
	if (token==ERROR) {
          printf("error in define_values\n");
	  exit(1);
        }
  }
/*printf("{END_INPUT}");*/
}

void
def_acc()
{
    def_named_variable("particles",T_DOUBLE_2);
    def_named_variable("energy",T_DOUBLE_2);
    def_named_variable("charge_sign",T_DOUBLE);
    def_named_variable("f_rep",T_DOUBLE);
    def_named_variable("n_b",T_INT);
    def_named_variable("sigma_x",T_DOUBLE_2);
    def_named_variable("sigma_y",T_DOUBLE_2);
    def_named_variable("sigma_z",T_DOUBLE_2);
    def_named_variable("emitt_x",T_DOUBLE_2);
    def_named_variable("emitt_y",T_DOUBLE_2);
    def_named_variable("beta_x",T_DOUBLE_2);
    def_named_variable("beta_y",T_DOUBLE_2);
    def_named_variable("offset_x",T_DOUBLE_MIRROR);
    def_named_variable("offset_y",T_DOUBLE_MIRROR);
    def_named_variable("offset_z",T_DOUBLE_MIRROR);
    def_named_variable("waist_x",T_DOUBLE_2);
    def_named_variable("waist_y",T_DOUBLE_2);
    def_named_variable("xi_x",T_DOUBLE_2);
    def_named_variable("xi_y",T_DOUBLE_2);
    def_named_variable("couple_xy",T_DOUBLE_2);
    def_named_variable("angle_phi",T_DOUBLE_2);
    def_named_variable("angle_x",T_DOUBLE_2);
    def_named_variable("angle_y",T_DOUBLE_2);
    def_named_variable("dist_x",T_INT_2);
    def_named_variable("dist_y",T_INT_2);
    def_named_variable("dist_z",T_INT_2);
    def_named_variable("trav_focus",T_INT_2);
    def_named_variable("espread",T_DOUBLE_2);
    def_named_variable("which_espread",T_INT_2);
/*scd 2.11.1998 DF*/
    def_named_variable("scale_step",T_DOUBLE_2);

#ifdef TWOBEAM
    def_named_variable("twobeam",T_INT);
    def_named_variable("particles_2",T_DOUBLE_2);
    def_named_variable("energy_2",T_DOUBLE_2);
    def_named_variable("charge_sign_2",T_DOUBLE);
    def_named_variable("sigma_x_2",T_DOUBLE_2);
    def_named_variable("sigma_y_2",T_DOUBLE_2);
    def_named_variable("sigma_z_2",T_DOUBLE_2);
    def_named_variable("emitt_x_2",T_DOUBLE_2);
    def_named_variable("emitt_y_2",T_DOUBLE_2);
    def_named_variable("beta_x_2",T_DOUBLE_2);
    def_named_variable("beta_y_2",T_DOUBLE_2);
    def_named_variable("offset_x_2",T_DOUBLE_MIRROR);
    def_named_variable("offset_y_2",T_DOUBLE_MIRROR);
    def_named_variable("offset_z_2",T_DOUBLE_MIRROR);
    def_named_variable("waist_x_2",T_DOUBLE_2);
    def_named_variable("waist_y_2",T_DOUBLE_2);
    def_named_variable("angle_phi_2",T_DOUBLE_2);
    def_named_variable("angle_x_2",T_DOUBLE_2);
    def_named_variable("angle_y_2",T_DOUBLE_2);
    def_named_variable("dist_x_2",T_INT_2);
    def_named_variable("dist_y_2",T_INT_2);
    def_named_variable("dist_z_2",T_INT_2);
    def_named_variable("trav_focus_2",T_INT_2);
    def_named_variable("espread_2",T_DOUBLE_2);
    def_named_variable("which_espread_2",T_INT_2);
#endif

    define_values("{xi_x=0.0;xi_y=0.0;}");
    define_values("{offset_x=0.0;offset_y=0.0;offset_z=0.0;}");
    define_values("{waist_x=0.0;waist_y=0.0;couple_xy=0.0;}");
    define_values("{angle_phi=0.0;angle_x=0.0;angle_y=0.0;}");
    define_values("{beta_x=-1.0;beta_y=-1.0;}");
    define_values("{sigma_x=-1.0;sigma_y=-1.0;}");
    define_values("{emitt_x=-1.0;emitt_y=-1.0;}");
    define_values("{charge_sign=-1;espread=1e-3;which_espread=1;}");
/*scd 2.11.1998 DF*/
    define_values("{scale_step=1.0;}");

#ifdef TWOBEAM
    define_values("{twobeam=0;}");
    define_values("{offset_x_2=0.0;offset_y_2=0.0;offset_z_2=0.0;}");
    define_values("{waist_x_2=0.0;waist_y_2=0.0;}");
    define_values("{angle_phi_2=0.0;angle_x_2=0.0;angle_y_2=0.0;}");
    define_values("{beta_x_2=-1.0;beta_y_2=-1.0;}");
    define_values("{sigma_x_2=-1.0;sigma_y_2=-1.0;}");
    define_values("{emitt_x_2=-1.0;emitt_y_2=-1.0;}");
    define_values("{charge_sign_2=-1;espread_2=1e-3;which_espread_2=1;}");
#endif
}

void
def_param()
{
  def_named_variable("silent",T_INT);
  def_named_variable("n_x",T_INT);
  def_named_variable("n_y",T_INT);
  def_named_variable("n_z",T_INT);
  def_named_variable("n_t",T_INT);
  def_named_variable("n_m",T_INT_2);
  def_named_variable("cut_x",T_DOUBLE);
  def_named_variable("cut_y",T_DOUBLE);
  def_named_variable("cut_z",T_DOUBLE);
  def_named_variable("integration_method",T_INT);
  def_named_variable("force_symmetric",T_INT);
  def_named_variable("charge_symmetric",T_INT);
  def_named_variable("do_size_log",T_INT);

  def_named_variable("rndm_load",T_INT);
  def_named_variable("rndm_save",T_INT);

  def_named_variable("do_photons",T_INT_2);
  def_named_variable("photon_ratio",T_DOUBLE);
  def_named_variable("ecm_min",T_DOUBLE);
  def_named_variable("ecm_min_gg",T_DOUBLE);

  def_named_variable("do_hadrons",T_INT);
  def_named_variable("store_hadrons",T_INT);
  def_named_variable("hadron_ratio",T_DOUBLE);

  def_named_variable("do_jets",T_INT);
  def_named_variable("jet_ptmin",T_DOUBLE);
  def_named_variable("jet_ratio",T_DOUBLE);
  def_named_variable("store_jets",T_INT);
  def_named_variable("jet_log",T_INT);

  def_named_variable("do_pairs",T_INT);
  def_named_variable("do_muons",T_INT);
  def_named_variable("do_coherent",T_INT);
  def_named_variable("track_pairs",T_INT);
  def_named_variable("store_pairs",T_INT);
  def_named_variable("pair_ecut",T_DOUBLE);
  def_named_variable("pair_step",T_DOUBLE);
  def_named_variable("beam_size",T_INT);
  def_named_variable("beam_size_scale",T_DOUBLE);
  def_named_variable("ext_field",T_INT);
  def_named_variable("pair_ratio",T_DOUBLE);
  def_named_variable("pair_q2",T_INT);
  def_named_variable("grids",T_INT);
  def_named_variable("beam_pair",T_INT);

  def_named_variable("electron_ratio",T_DOUBLE);
  def_named_variable("do_eloss",T_INT);
  def_named_variable("do_espread",T_INT);
  def_named_variable("do_isr",T_INT);
  def_named_variable("do_lumi",T_INT);

  def_named_variable("do_bhabha",T_INT);
  def_named_variable("bhabha_scal",T_DOUBLE);

  def_named_variable("do_prod",T_INT);
  def_named_variable("prod_e",T_DOUBLE);
  def_named_variable("prod_scal",T_DOUBLE);

  def_named_variable("load_events",T_INT);

  def_named_variable("do_lumi_ee_2",T_INT);
  def_named_variable("hist_ee_bins",T_INT);
  def_named_variable("hist_ee_min",T_DOUBLE);
  def_named_variable("hist_ee_max",T_DOUBLE);
  def_named_variable("hist_espec_bins",T_INT);
  def_named_variable("hist_espec_min",T_DOUBLE);
  def_named_variable("hist_espec_max",T_DOUBLE);
  def_named_variable("lumi_ee_2_n",T_INT);
  def_named_variable("lumi_ee_2_min",T_DOUBLE);
  def_named_variable("lumi_ee_2_max",T_DOUBLE);
  def_named_variable("do_lumi_eg_2",T_INT);
  def_named_variable("lumi_eg_2_n",T_INT);
  def_named_variable("lumi_eg_2_min",T_DOUBLE);
  def_named_variable("lumi_eg_2_max",T_DOUBLE);
  def_named_variable("do_lumi_ge_2",T_INT);
  def_named_variable("lumi_ge_2_n",T_INT);
  def_named_variable("lumi_ge_2_min",T_DOUBLE);
  def_named_variable("lumi_ge_2_max",T_DOUBLE);
  def_named_variable("do_lumi_gg_2",T_INT);
  def_named_variable("lumi_gg_2_n",T_INT);
  def_named_variable("lumi_gg_2_min",T_DOUBLE);
  def_named_variable("lumi_gg_2_max",T_DOUBLE);

  def_named_variable("beam_vx_min",T_DOUBLE);
  def_named_variable("beam_vx_max",T_DOUBLE);
  def_named_variable("beam_vx_interval",T_INT);
  def_named_variable("beam_vy_min",T_DOUBLE);
  def_named_variable("beam_vy_max",T_DOUBLE);
  def_named_variable("beam_vy_interval",T_INT);

  def_named_variable("do_cross",T_INT);
  def_named_variable("do_compt",T_INT);
  def_named_variable("do_compt_phot",T_INT);
  def_named_variable("compt_scale",T_DOUBLE);
  def_named_variable("compt_x_min",T_DOUBLE);
  def_named_variable("compt_emax",T_DOUBLE);
  def_named_variable("num_lumi",T_INT);
  def_named_variable("num_lumi_eg",T_INT);
  def_named_variable("num_lumi_gg",T_INT);
  def_named_variable("lumi_p",T_DOUBLE);
  def_named_variable("lumi_p_eg",T_DOUBLE);
  def_named_variable("lumi_p_gg",T_DOUBLE);
  def_named_variable("store_beam",T_INT);
  def_named_variable("load_beam",T_INT);
  def_named_variable("store_photons",T_INT);
  def_named_variable("load_photons",T_INT);
  def_named_variable("jet_pythia",T_INT);

  def_named_variable("do_dump",T_INT);
  def_named_variable("dump_step",T_INT);
  def_named_variable("dump_particle",T_INT);

#ifdef TWOBEAM
  def_named_variable("n2_m",T_INT_2);
#endif

  define_values("{silent=1;integration_method=2;rndm_save=1;rndm_load=1;}");
  define_values("{do_eloss=1;do_isr=0;do_espread=0;do_cross=0;}");
  define_values("{do_prod=0;prod_e=0.0;prod_scal=1e-29;load_events=0;}");
  define_values("{load_beam=0;store_beam=0;load_photons=0;store_photons=0;}");
  define_values("{do_photons=0;photon_ratio=1.0;electron_ratio=1.0;}");
  define_values("{do_hadrons=0;hadron_ratio=1e5;}");
  define_values("{do_jets=0;jet_ratio=1e5;store_jets=0;jet_log=1;}");
  define_values("{jet_ptmin=2.0;}");
  define_values("{beam_pair=0;n_b=1;f_rep=1;}");
  define_values("{do_pairs=0;track_pairs=0;pair_ratio=1.0;pair_q2=2;}");
  define_values("{pair_ecut=5e-3;pair_step=1.0;do_muons=0;}");
  define_values("{beam_size=1;beam_size_scale=1.0;ext_field=0;store_pairs=1;}");
  define_values("{do_compt=0;do_compt_phot=0;compt_x_min=1.0;}");
  define_values("{compt_emax=100.0;do_coherent=0;compt_scale=1.0;}");
  define_values("{hist_ee_bins=200;hist_ee_min=0.0;hist_ee_max=-1.0;}");
  define_values("{hist_espec_bins=-1;hist_espec_min=0.0;}");
  define_values("{hist_espec_max=500.001;}");
  define_values("{do_lumi_ee_2=0;num_lumi=0;do_size_log=0;}");
  define_values("{lumi_ee_2_n=100;lumi_ee_2_min=0.0;lumi_ee_2_max=1.0;}");
  define_values("{do_lumi_eg_2=0;}");
  define_values("{lumi_eg_2_n=100;lumi_eg_2_min=0.0;lumi_eg_2_max=1.0;}");
  define_values("{do_lumi_ge_2=0;}");
  define_values("{lumi_ge_2_n=100;lumi_ge_2_min=0.0;lumi_ge_2_max=1.0;}");
  define_values("{do_lumi_gg_2=0;}");
  define_values("{do_lumi=0;num_lumi=100000;lumi_p=1e-29;}");
  define_values("{lumi_gg_2_n=100;lumi_gg_2_min=0.0;lumi_gg_2_max=1.0;}");
  define_values("{do_bhabha=0;bhabha_scal=0.0;ecm_min=0.0;ecm_min_gg=0.0;}");
  define_values("{cut_x=-1.0;cut_y=-1.0;cut_z=-1.0;}");
  define_values("{jet_pythia=0;}");
  define_values("{beam_vx_min=-1e-3;beam_vx_max=1e-3;beam_vx_interval=200;}");
  define_values("{beam_vy_min=-1e-3;beam_vy_max=1e-3;beam_vy_interval=200;}");
  define_values("{do_dump=0;dump_step=1;dump_particle=1;}");
}

void
input_values(char *buffer)
{
  TOKEN token;
  VALUE value;
  store_buffer(buffer);
  if ((token=read_token(&value))!=START){
	printf("error: no proper start\n");
	exit(1);
  }
  while((token=read_command())!=END_INPUT) {
/*	printf("{");
	print_token(token);
	printf("}");*/
	if (token==ERROR) {
	  printf("error while scanning input list\n");
	  exit(1);
        }
  }
}

struct struct_datei
{
  FILE *file;
  char buffer[max_buffer];
  char *point,*line_end;
};

typedef struct struct_datei DATEI;

int file_open(DATEI *datei,char *name,char *type)
{
  if (name) {
    if(!(datei->file=fopen(name,type))) return 0;
    datei->buffer[0]='\0';
    datei->point=datei->buffer;
    datei->line_end=datei->buffer;
    return 1;
  }
  else {
    datei->file=stdin;
    datei->buffer[0]='\0';
    datei->point=datei->buffer;
    datei->line_end=datei->buffer;
    return 1;
  }
}

void file_close(DATEI *datei)
{
  if (datei->file==stdin) {
  }
  else {
    fclose(datei->file);
  }
}

void file_rewind(DATEI *datei)
{
  rewind(datei->file);
}

int file_read_line(DATEI *datei)
{
  if(!(datei->point=fgets(datei->buffer,max_buffer,datei->file))) return 0;
  datei->line_end=datei->point+strlen(datei->buffer);
  *(--(datei->line_end))='\0';
  return 1;
}

int file_find_word_line(DATEI *datei,char *word)
{
  if(!(datei->point=strstr(datei->point,word))) return 0;
  datei->point+=strlen(word);
  return 1;
}

int file_find_word(DATEI *datei,char *word)
{
  if (file_find_word_line(datei,word)) return 1;
  while(file_read_line(datei))
    {
      if(file_find_word_line(datei,word)) return 1;
    }
  return 0;
}

int file_skip_space_line(DATEI *datei)
{
  if(datei->point==NULL) return 0;
  (datei->point)--;
  while((++(datei->point)) < datei->line_end)
    if ((*(datei->point))!=' ') return 1;
  return 0;
}

int file_skip_space(DATEI *datei)
{
  if (file_skip_space_line(datei)) return 1;
  while(file_read_line(datei))
    if (file_skip_space_line(datei)) return 1;
 return 0;
}

int file_next_word(DATEI *datei,char *word)
{
  file_skip_space(datei);
  if (strncmp(datei->point,word,strlen(word))!=0) return 0;
  datei->point+=strlen(word);
  return 1;
}

double file_get_double(DATEI *datei)
{
  if(!file_skip_space(datei)) return 0;
  return strtod(datei->point,&(datei->point));
}

int file_read_double(DATEI *datei,double *x)
{
  if(!file_skip_space(datei)) return 0;
  *x=strtod(datei->point,&(datei->point));
  return 1;
}

int file_read_int(DATEI *datei,int *n)
{
  if(!file_skip_space(datei)) return 0;
  *n=strtol(datei->point,&(datei->point),10);
  return 1;
}

#define stpcpy(a,b) (strcpy(a,b)+strlen(b))

int file_read_until(DATEI *datei,char *end,char *buff,int n_max)
{
  char *tmp;
  int len;
  buff[0]='\0';
  len=strlen(end);
  while((tmp=strstr(datei->point,end))==NULL)
    {
      if((n_max-=strlen(datei->point))>1)
	{
	  buff=stpcpy(buff,datei->point);
	}
      else
	{
	  return 0;
	}
      file_read_line(datei);
    }
  if((n_max-=(tmp-datei->point)+len)>1)
      {
	buff=strncat(buff,datei->point,tmp-datei->point+len);
      }
  else{
    }
  datei->point=tmp+len;
  return 1;
}

int file_read_braces(DATEI *datei,char *begin,char *end,char *buff,int n_max)
{
  if(file_skip_space(datei)==0) return 0;
  if(strncmp(datei->point,begin,strlen(begin))!=0) return 0;
  file_read_until(datei,end,buff,n_max);
  return 1;
}

int file_find_braces(DATEI *datei,char *begin,char *end,char *buff,int n_max)
{
  if (!file_find_word(datei,begin)) return 0;
  if (!file_read_until(datei,end,buff,n_max)) return 0;
  return 1;
}

//#include "hist.d"
#include "histogram.h"

int string_get_int(char *line,char *word,int *n)
{
  char *pointer;

  if ((pointer=strstr(line,word))==NULL) return 0;
  pointer+=strlen(word);
  while (*pointer==' ') pointer++;
  if (*pointer!='=') return 0;
  pointer++;
  while (*pointer==' ') pointer++;
  *n=strtol(pointer,&pointer,10);
  return 1;
}

int string_get_float(char *line,char *word,float *x)
{
  char *pointer;

  if ((pointer=strstr(line,word))==NULL) return 0;
  pointer+=strlen(word);
  while (*pointer==' ') pointer++;
  if (*pointer!='=') return 0;
  pointer++;
  while (*pointer==' ') pointer++;
  *x=strtod(pointer,&pointer);
  return 1;
}

int string_get_list(char *line,char *word,float *x,int *n)
{
  char *pointer;
  int i=0;

  if ((pointer=strstr(line,word))==NULL) return 0;
  pointer+=strlen(word);
  while (*pointer==' ') pointer++;
  if (*pointer!='=') return 0;
  pointer++;
  while (*pointer==' ') pointer++;
  if(*pointer!='(') return 0;
  pointer++;
  while (*pointer!=')'){
     x[i]=strtod(pointer,&pointer);
     i++;
     pointer++;
  }
  *n=i;
  return 1;
}

int file_read_histogram(DATEI *datei,char *name,HISTOGRAM *histogram)
{
  int i,n;
  double x,x2;
  do
    {
      if (file_find_word(datei,"$HISTOGRAM::")==0) return 0;
    }
  while(strstr(datei->point,name)!=datei->point);
  file_read_line(datei);
  string_get_int(datei->point,"type",&histogram->type);
  string_get_int(datei->point,"nc",&histogram->nc);
  string_get_float(datei->point,"xmin",&histogram->xmin);
  string_get_float(datei->point,"xmax",&histogram->xmax);
  if (histogram->type==2)
    {
      histogram->xmin=exp(histogram->xmin);
      histogram->xmax=exp(histogram->xmax);
    }
  histogram_make(histogram,histogram->type,histogram->xmin,histogram->xmax,
		 histogram->nc,name);
  file_read_line(datei);
  switch(histogram->type)
    {
    case 1:
    case 2:
      for (i=0;i<histogram->nc;i++)
	{
	  file_read_double(datei,&x);
	  histogram_set_channel_y(histogram,i,x);
	}
      break;
    case 3:
      for (i=0;i<histogram->nc;i++)
	{
	  file_read_double(datei,&x);
	  file_read_int(datei,&n);
	  histogram_set_channel_y(histogram,i,x);
	  histogram_set_channel_count(histogram,i,n);
	}
      break;
    case 4:
      for (i=0;i<histogram->nc;i++)
	{
	  file_read_int(datei,&n);
	  file_read_double(datei,&x);
	  file_read_double(datei,&x2);
	  histogram_set_channel_count(histogram,i,n);
	  histogram_set_channel_y(histogram,i,x);
	  histogram_set_channel_y2(histogram,i,x2);
	}
      break;
    }
  return 1;
}

void file_list_histograms(DATEI *datei)
{
  while (file_find_word(datei,"$HISTOGRAM::")!=0)
    printf("%s\n",datei->point);
}

/*
int file_read_next_word(DATEI *datei,char *word)
{
}
*/

#ifdef MUON_BEAM

void acc_test(float ebeam,float *emitt,float *beta,float *sigma)
{
/*  emitt*=0.511e-3/ebeam;*/
printf("acc: %g %g %g\n",*emitt,*beta,*sigma);
  if (*beta>0.0)
    {
	if (*emitt>0.0){
           *sigma=sqrt(*emitt *MUMASS/ebeam * *beta)*1e9;
        }
	else{
	   *emitt=*sigma * *sigma *1e-18 / *beta * ebeam/MUMASS;
	}
    }
  else
    {
      *beta=(*sigma * *sigma*1e-18)/(*emitt*MUMASS/ebeam);
    }
}

#else

void acc_test(float ebeam,float *emitt,float *beta,float *sigma)
{
/*  emitt*=0.511e-3/ebeam;*/
printf("acc: %g %g %g\n",*emitt,*beta,*sigma);
  if (*beta>0.0)
    {
	if (*emitt>0.0){
           *sigma=sqrt(*emitt *EMASS/ebeam * *beta)*1e9;
        }
	else{
	   *emitt=*sigma * *sigma *1e-18 / *beta * ebeam/EMASS;
	}
    }
  else
    {
      *beta=(*sigma * *sigma*1e-18)/(*emitt*EMASS/ebeam);
    }
}

#endif

#ifdef TWOBEAM
int file_read_accelerator_i(BEAM_PARAMETERS *beam1,BEAM_PARAMETERS *beam2,
	                    BEAM_PARAMETERS *beam1_2,BEAM_PARAMETERS *beam2_2,
                            GRID *grid,SWITCHES *switches)
#else
int file_read_accelerator_i(BEAM_PARAMETERS *beam1,BEAM_PARAMETERS *beam2,
	GRID *grid,SWITCHES *switches)
#endif
{
  VALUE value;
  float x;
  int m;

  get_named_variable("energy.1",&value);
  beam1->ebeam=CONTENTS(value);
  get_named_variable("energy.2",&value);
  beam2->ebeam=CONTENTS(value);

  get_named_variable("particles.1",&value);
  beam1->n_particles=CONTENTS(value)*1e10;
  get_named_variable("particles.2",&value);
  beam2->n_particles=CONTENTS(value)*1e10;

  get_named_variable("emitt_x.1",&value);
  beam1->em_x=CONTENTS(value)*1e-6;
  get_named_variable("emitt_x.2",&value);
  beam2->em_x=CONTENTS(value)*1e-6;

  get_named_variable("emitt_y.1",&value);
  beam1->em_y=CONTENTS(value)*1e-6;
  get_named_variable("emitt_y.2",&value);
  beam2->em_y=CONTENTS(value)*1e-6;

  get_named_variable("beta_x.1",&value);
  beam1->beta_x=CONTENTS(value)*1e-3;
  get_named_variable("beta_x.2",&value);
  beam2->beta_x=CONTENTS(value)*1e-3;

  get_named_variable("beta_y.1",&value);
  beam1->beta_y=CONTENTS(value)*1e-3;
  get_named_variable("beta_y.2",&value);
  beam2->beta_y=CONTENTS(value)*1e-3;

  get_named_variable("sigma_x.1",&value);
  beam1->sigma_x=CONTENTS(value);
  get_named_variable("sigma_x.2",&value);
  beam2->sigma_x=CONTENTS(value);

  get_named_variable("sigma_y.1",&value);
  beam1->sigma_y=CONTENTS(value);
  get_named_variable("sigma_y.2",&value);
  beam2->sigma_y=CONTENTS(value);

  get_named_variable("sigma_z.1",&value);
  beam1->sigma_z=CONTENTS(value)*1e3;
  get_named_variable("sigma_z.2",&value);
  beam2->sigma_z=CONTENTS(value)*1e3;

  get_named_variable("dist_z.1",&value);
  beam1->dist_z=CONTENTS(value);
  get_named_variable("dist_z.2",&value);
  beam2->dist_z=CONTENTS(value);

  get_named_variable("dist_x.1",&value);
  beam1->dist_x=CONTENTS(value);
  get_named_variable("dist_x.2",&value);
  beam2->dist_x=CONTENTS(value);

  get_named_variable("dist_y.1",&value);
  beam1->dist_y=CONTENTS(value);
  get_named_variable("dist_y.2",&value);
  beam2->dist_y=CONTENTS(value);

  get_named_variable("trav_focus.1",&value);
  beam1->trav_focus=CONTENTS(value);
  get_named_variable("trav_focus.2",&value);
  beam2->trav_focus=CONTENTS(value);

  get_named_variable("offset_x.1",&value);
  beam1->offset_x=CONTENTS(value);
  get_named_variable("offset_x.2",&value);
  beam2->offset_x=CONTENTS(value);

  get_named_variable("offset_y.1",&value);
  beam1->offset_y=CONTENTS(value);
  get_named_variable("offset_y.2",&value);
  beam2->offset_y=CONTENTS(value);

  get_named_variable("offset_z.1",&value);
  beam1->offset_z=CONTENTS(value)*1e3;
  get_named_variable("offset_z.2",&value);
  beam2->offset_z=CONTENTS(value)*1e3;

  get_named_variable("waist_x.1",&value);
  beam1->waist_x=CONTENTS(value)*1e3;
  get_named_variable("waist_x.2",&value);
  beam2->waist_x=CONTENTS(value)*1e3;

  get_named_variable("waist_y.1",&value);
  beam1->waist_y=CONTENTS(value)*1e3;
  get_named_variable("waist_y.2",&value);
  beam2->waist_y=CONTENTS(value)*1e3;

  get_named_variable("couple_xy.1",&value);
  beam1->couple_xy=CONTENTS(value)*1e3;
  get_named_variable("couple_xy.2",&value);
  beam2->couple_xy=CONTENTS(value)*1e3;

  get_named_variable("xi_x.1",&value);
  beam1->xi_x=CONTENTS(value)*1e3;
  get_named_variable("xi_x.2",&value);
  beam2->xi_x=CONTENTS(value)*1e3;

  get_named_variable("xi_y.1",&value);
  beam1->xi_y=CONTENTS(value);
  get_named_variable("xi_y.2",&value);
  beam2->xi_y=CONTENTS(value);

  get_named_variable("angle_phi.1",&value);
  beam1->phi_angle=CONTENTS(value);
  get_named_variable("angle_phi.2",&value);
  beam2->phi_angle=CONTENTS(value);

  get_named_variable("angle_x.1",&value);
  beam1->x_angle=CONTENTS(value);
  get_named_variable("angle_x.2",&value);
  beam2->x_angle=CONTENTS(value);

  get_named_variable("angle_y.1",&value);
  beam1->y_angle=CONTENTS(value);
  get_named_variable("angle_y.2",&value);
  beam2->y_angle=CONTENTS(value);

  get_named_variable("n_b",&value);
  beam1->bunches_per_train=CONTENTS(value);
  beam2->bunches_per_train=CONTENTS(value);
  switches->n_b=CONTENTS(value);

  get_named_variable("f_rep",&value);
  beam1->frep=CONTENTS(value);
  beam2->frep=CONTENTS(value);
  switches->f_rep=CONTENTS(value);

  get_named_variable("charge_sign",&value);
  switches->charge_sign=CONTENTS(value);
  switches->charge_sign_0=0.0;
  if(switches->charge_sign>0.0){
    switches->charge_sign_0=1.0;
  }
  if(switches->charge_sign<0.0){
    switches->charge_sign_0=-1.0;
  }

  acc_test(beam1->ebeam,&(beam1->em_x),&(beam1->beta_x),&(beam1->sigma_x));
  double_to_value(beam1->em_x*1e6,&value);
  set_named_variable("emitt_x.1",value);
  double_to_value(beam1->beta_x*1e3,&value);
  set_named_variable("beta_x.1",value);
  double_to_value(beam1->sigma_x,&value);
  set_named_variable("sigma_x.1",value);

  acc_test(beam2->ebeam,&(beam2->em_x),&(beam2->beta_x),&(beam2->sigma_x));
  double_to_value(beam2->em_x*1e6,&value);
  set_named_variable("emitt_x.2",value);
  double_to_value(beam2->beta_x*1e3,&value);
  set_named_variable("beta_x.2",value);
  double_to_value(beam2->sigma_x,&value);
  set_named_variable("sigma_x.2",value);

  acc_test(beam1->ebeam,&(beam1->em_y),&(beam1->beta_y),&(beam1->sigma_y));
  double_to_value(beam1->em_y*1e6,&value);
  set_named_variable("emitt_y.1",value);
  double_to_value(beam1->beta_y*1e3,&value);
  set_named_variable("beta_y.1",value);
  double_to_value(beam1->sigma_y,&value);
  set_named_variable("sigma_y.1",value);

  acc_test(beam2->ebeam,&(beam2->em_y),&(beam2->beta_y),&(beam2->sigma_y));
  double_to_value(beam2->em_y*1e6,&value);
  set_named_variable("emitt_y.2",value);
  double_to_value(beam2->beta_y*1e3,&value);
  set_named_variable("beta_y.2",value);
  double_to_value(beam2->sigma_y,&value);
  set_named_variable("sigma_y.2",value);

#ifdef TWOBEAM
  get_named_variable("twobeam",&value);
  switches->twobeam=CONTENTS(value);
  if (switches->twobeam){
      get_named_variable("energy_2.1",&value);
      beam1_2->ebeam=CONTENTS(value);
      get_named_variable("energy_2.2",&value);
      beam2_2->ebeam=CONTENTS(value);

      get_named_variable("particles_2.1",&value);
      beam1_2->n_particles=CONTENTS(value)*1e10;
      get_named_variable("particles_2.2",&value);
      beam2_2->n_particles=CONTENTS(value)*1e10;

      get_named_variable("emitt_x_2.1",&value);
      beam1_2->em_x=CONTENTS(value)*1e-6;
      get_named_variable("emitt_x_2.2",&value);
      beam2_2->em_x=CONTENTS(value)*1e-6;

      get_named_variable("emitt_y_2.1",&value);
      beam1_2->em_y=CONTENTS(value)*1e-6;
      get_named_variable("emitt_y_2.2",&value);
      beam2_2->em_y=CONTENTS(value)*1e-6;

      get_named_variable("beta_x_2.1",&value);
      beam1_2->beta_x=CONTENTS(value)*1e-3;
      get_named_variable("beta_x_2.2",&value);
      beam2_2->beta_x=CONTENTS(value)*1e-3;

      get_named_variable("beta_y_2.1",&value);
      beam1_2->beta_y=CONTENTS(value)*1e-3;
      get_named_variable("beta_y_2.2",&value);
      beam2_2->beta_y=CONTENTS(value)*1e-3;

      get_named_variable("sigma_x_2.1",&value);
      beam1_2->sigma_x=CONTENTS(value);
      get_named_variable("sigma_x_2.2",&value);
      beam2_2->sigma_x=CONTENTS(value);

      get_named_variable("sigma_y_2.1",&value);
      beam1_2->sigma_y=CONTENTS(value);
      get_named_variable("sigma_y_2.2",&value);
      beam2_2->sigma_y=CONTENTS(value);

      get_named_variable("sigma_z_2.1",&value);
      beam1_2->sigma_z=CONTENTS(value)*1e3;
      get_named_variable("sigma_z_2.2",&value);
      beam2_2->sigma_z=CONTENTS(value)*1e3;

      get_named_variable("dist_z_2.1",&value);
      beam1_2->dist_z=CONTENTS(value);
      get_named_variable("dist_z_2.2",&value);
      beam2_2->dist_z=CONTENTS(value);

      get_named_variable("dist_x_2.1",&value);
      beam1_2->dist_x=CONTENTS(value);
      get_named_variable("dist_x_2.2",&value);
      beam2_2->dist_x=CONTENTS(value);

      get_named_variable("dist_y_2.1",&value);
      beam1_2->dist_y=CONTENTS(value);
      get_named_variable("dist_y_2.2",&value);
      beam2_2->dist_y=CONTENTS(value);

      get_named_variable("trav_focus_2.1",&value);
      beam1_2->trav_focus=CONTENTS(value);
      get_named_variable("trav_focus_2.2",&value);
      beam2_2->trav_focus=CONTENTS(value);

      get_named_variable("offset_x_2.1",&value);
      beam1_2->offset_x=CONTENTS(value);
      get_named_variable("offset_x_2.2",&value);
      beam2_2->offset_x=CONTENTS(value);

      get_named_variable("offset_y_2.1",&value);
      beam1_2->offset_y=CONTENTS(value);
      get_named_variable("offset_y_2.2",&value);
      beam2_2->offset_y=CONTENTS(value);

      get_named_variable("offset_z_2.1",&value);
      beam1_2->offset_z=CONTENTS(value)*1e3;
      get_named_variable("offset_z_2.2",&value);
      beam2_2->offset_z=CONTENTS(value)*1e3;

      get_named_variable("waist_x_2.1",&value);
      beam1_2->waist_x=CONTENTS(value)*1e3;
      get_named_variable("waist_x_2.2",&value);
      beam2_2->waist_x=CONTENTS(value)*1e3;

      get_named_variable("waist_y_2.1",&value);
      beam1_2->waist_y=CONTENTS(value)*1e3;
      get_named_variable("waist_y_2.2",&value);
      beam2_2->waist_y=CONTENTS(value)*1e3;

      get_named_variable("angle_phi_2.1",&value);
      beam1_2->phi_angle=CONTENTS(value);
      get_named_variable("angle_phi_2.2",&value);
      beam2_2->phi_angle=CONTENTS(value);

      get_named_variable("angle_x_2.1",&value);
      beam1_2->x_angle=CONTENTS(value);
      get_named_variable("angle_x_2.2",&value);
      beam2_2->x_angle=CONTENTS(value);

      get_named_variable("angle_y_2.1",&value);
      beam1_2->y_angle=CONTENTS(value);
      get_named_variable("angle_y_2.2",&value);
      beam2_2->y_angle=CONTENTS(value);

      get_named_variable("charge_sign_2",&value);
      switches->charge_sign_2=CONTENTS(value);

      acc_test(beam1_2->ebeam,&(beam1_2->em_x),&(beam1_2->beta_x),
               &(beam1_2->sigma_x));
      double_to_value(beam1_2->em_x*1e6,&value);
      set_named_variable("emitt_x_2.1",value);
      double_to_value(beam1_2->beta_x*1e3,&value);
      set_named_variable("beta_x_2.1",value);
      double_to_value(beam1_2->sigma_x,&value);
      set_named_variable("sigma_x_2.1",value);

      acc_test(beam2_2->ebeam,&(beam2_2->em_x),&(beam2_2->beta_x),
               &(beam2_2->sigma_x));
      double_to_value(beam2_2->em_x*1e6,&value);
      set_named_variable("emitt_x.2",value);
      double_to_value(beam2_2->beta_x*1e3,&value);
      set_named_variable("beta_x.2",value);
      double_to_value(beam2_2->sigma_x,&value);
      set_named_variable("sigma_x.2",value);

      acc_test(beam1_2->ebeam,&(beam1_2->em_y),&(beam1_2->beta_y),
               &(beam1_2->sigma_y));
      double_to_value(beam1_2->em_y*1e6,&value);
      set_named_variable("emitt_y_2.1",value);
      double_to_value(beam1_2->beta_y*1e3,&value);
      set_named_variable("beta_y_2.1",value);
      double_to_value(beam1_2->sigma_y,&value);
      set_named_variable("sigma_y_2.1",value);

      acc_test(beam2_2->ebeam,&(beam2_2->em_y),&(beam2_2->beta_y),
               &(beam2_2->sigma_y));
      double_to_value(beam2_2->em_y*1e6,&value);
      set_named_variable("emitt_y_2.2",value);
      double_to_value(beam2_2->beta_y*1e3,&value);
      set_named_variable("beta_y_2.2",value);
      double_to_value(beam2_2->sigma_y,&value);
      set_named_variable("sigma_y_2.2",value);
  }
#endif
}

#ifdef TWOBEAM
int file_read_accelerator(DATEI *datei,char *name,BEAM_PARAMETERS *beam1,
			  BEAM_PARAMETERS *beam2,BEAM_PARAMETERS *beam1_2,
                          BEAM_PARAMETERS *beam2_2,GRID *grid,
                          SWITCHES *switches)
#else
int file_read_accelerator(DATEI *datei,char *name,BEAM_PARAMETERS *beam1,
			  BEAM_PARAMETERS *beam2,GRID *grid,SWITCHES *switches)
#endif
{
  static int n=max_buffer2;
  char buffer[max_buffer2];
  def_acc();
  do
    {
      do
	{
	  if (file_find_word(datei,"$ACCELERATOR::")==0) return 0;
	}
      while(file_next_word(datei,name)==0);
    }
  while(!file_read_braces(datei,"{","}",buffer,n)) ;

  input_values(buffer);
#ifdef TWOBEAM
  file_read_accelerator_i(beam1,beam2,beam1_2,beam2_2,grid,switches);
#else
  file_read_accelerator_i(beam1,beam2,grid,switches);
#endif
  return 1;
}

/*
#ifdef TWOBEAM
int file_read_accelerator2(DATEI *datei,BEAM_PARAMETERS *beam1,
			  BEAM_PARAMETERS *beam2,BEAM_PARAMETERS *beam1_2,
                          BEAM_PARAMETERS *beam2_2,GRID *grid,
                          SWITCHES *switches)
#else
int file_read_accelerator2(DATEI *datei,BEAM_PARAMETERS *beam1,
			  BEAM_PARAMETERS *beam2,GRID *grid,SWITCHES *switches)
#endif
{
  char buffer[10000];
  file_read_braces(datei,"{","}",buffer,10000);
#ifdef TWOBEAM
  file_read_accelerator_i(beam1,beam2,beam1_2,beam2_2,grid,switches);
#else
  file_read_accelerator_i(beam1,beam2,grid,switches);
#endif
  return 1;
}
*/

#ifdef TWOBEAM
int read_file_parameters_i(BEAM_PARAMETERS *beam1,BEAM_PARAMETERS *beam2,
                           BEAM_PARAMETERS *beam1_2,BEAM_PARAMETERS *beam2_2,
                           GRID *grid,SWITCHES *switches)
#else
int read_file_parameters_i(BEAM_PARAMETERS *beam1,
	                   BEAM_PARAMETERS *beam2,GRID *grid,
                           SWITCHES *switches)
#endif
{
//  const int n=max_buffer2;
//  char buffer[max_buffer2];
  float tmp;
  VALUE value;

  get_named_variable("integration_method",&value);
  switches->integration_method=CONTENTS(value);

  get_named_variable("silent",&value);
  switches->silent=CONTENTS(value);

  get_named_variable("n_x",&value);
  grid->n_cell_x=CONTENTS(value);

  get_named_variable("n_y",&value);
  grid->n_cell_y=CONTENTS(value);

  get_named_variable("n_z",&value);
  grid->n_cell_z=CONTENTS(value);

  get_named_variable("n_t",&value);
  grid->timestep=CONTENTS(value);

  get_named_variable("n_m.1",&value);
  grid->n_m_1=CONTENTS(value);
  beam1->n_macro=CONTENTS(value);
  get_named_variable("n_m.2",&value);
  grid->n_m_2=CONTENTS(value);
  beam2->n_macro=CONTENTS(value);

#ifdef TWOBEAM
  get_named_variable("n2_m.1",&value);
  beam1_2->n_macro=CONTENTS(value);
  get_named_variable("n2_m.2",&value);
  beam2_2->n_macro=CONTENTS(value);
#endif

  get_named_variable("cut_x",&value);
  if(CONTENTS(value)>0.0){
      grid->cut_x=CONTENTS(value);
  }
  else{
      get_named_variable("sigma_x.1",&value);
      tmp=CONTENTS(value);
      get_named_variable("sigma_x.2",&value);
      if(tmp<CONTENTS(value)) tmp=CONTENTS(value);
      grid->cut_x=3.0*tmp;
  }

  get_named_variable("cut_y",&value);
  if(CONTENTS(value)>0.0){
      grid->cut_y=CONTENTS(value);
  }
  else{
      get_named_variable("sigma_y.1",&value);
      tmp=CONTENTS(value);
      get_named_variable("sigma_y.2",&value);
      if(tmp<CONTENTS(value)) tmp=CONTENTS(value);
      grid->cut_y=3.0*tmp;
  }

  get_named_variable("cut_z",&value);
  if(CONTENTS(value)>0.0){
      grid->cut_z=CONTENTS(value)*1e3;
  }
  else{
      get_named_variable("sigma_z.1",&value);
      tmp=CONTENTS(value);
      get_named_variable("sigma_z.2",&value);
      if(tmp<CONTENTS(value)) tmp=CONTENTS(value);
      grid->cut_z=3.0*tmp*1e3;
  }

/*scd 2.11.1998 DF*/
  get_named_variable("scale_step.1",&value);
  grid->scal_step[1]=CONTENTS(value);
  get_named_variable("scale_step.2",&value);
  grid->scal_step[0]=CONTENTS(value);

  get_named_variable("grids",&value);
  switches->extra_grids=CONTENTS(value)-1;
  if (switches->extra_grids<0) switches->extra_grids=0;

  get_named_variable("load_photons",&value);
  switches->load_photon=CONTENTS(value);

  get_named_variable("store_photons",&value);
  switches->write_photons=CONTENTS(value);

  get_named_variable("do_photons.1",&value);
  switches->store_photons[1]=CONTENTS(value);
  get_named_variable("do_photons.2",&value);
  switches->store_photons[2]=CONTENTS(value);
  
  get_named_variable("do_eloss",&value);
  switches->do_eloss=CONTENTS(value);

  get_named_variable("ecm_min",&value);
  switches->ecm_min=CONTENTS(value);

  get_named_variable("ecm_min_gg",&value);
  switches->gg_smin=4*(CONTENTS(value))*(CONTENTS(value));

  get_named_variable("do_hadrons",&value);
  switches->do_hadrons=CONTENTS(value);

  get_named_variable("store_hadrons",&value);
  switches->store_hadrons=CONTENTS(value);

  get_named_variable("hadron_ratio",&value);
  switches->hadron_ratio=CONTENTS(value);

  get_named_variable("do_jets",&value);
  switches->do_jets=CONTENTS(value);

  get_named_variable("do_pairs",&value);
  switches->do_pairs=CONTENTS(value);

  get_named_variable("beam_pair",&value);
  switches->beam_pair=CONTENTS(value);

  get_named_variable("store_pairs",&value);
  switches->store_pairs=CONTENTS(value);

  get_named_variable("do_muons",&value);
  switches->do_muons=CONTENTS(value);

  get_named_variable("do_coherent",&value);
  switches->do_coherent=CONTENTS(value);

  get_named_variable("track_pairs",&value);
  switches->track_pairs=CONTENTS(value);

  get_named_variable("pair_ecut",&value);
  switches->pair_ecut=CONTENTS(value);

  get_named_variable("pair_step",&value);
  switches->pair_step=CONTENTS(value);

  get_named_variable("electron_ratio",&value);
  switches->electron_ratio=CONTENTS(value);

  get_named_variable("do_lumi",&value);
  switches->do_lumi=CONTENTS(value);

  get_named_variable("bhabha_scal",&value);
  switches->bhabha_scal=CONTENTS(value);

  get_named_variable("rndm_save",&value);
  switches->rndm_save=CONTENTS(value);

  get_named_variable("rndm_load",&value);
  switches->rndm_load=CONTENTS(value);

  get_named_variable("hist_ee_min",&value);
  switches->hist_ee_min=CONTENTS(value);

  get_named_variable("hist_ee_max",&value);
  switches->hist_ee_max=CONTENTS(value);

  get_named_variable("hist_ee_bins",&value);
  switches->hist_ee_bins=CONTENTS(value);

  get_named_variable("hist_espec_min",&value);
  switches->hist_espec_min=CONTENTS(value);

  get_named_variable("hist_espec_max",&value);
  switches->hist_espec_max=CONTENTS(value);

  get_named_variable("hist_espec_bins",&value);
  switches->hist_espec_bins=CONTENTS(value);

  get_named_variable("do_lumi_ee_2",&value);
  switches->do_lumi_ee_2=CONTENTS(value);

  get_named_variable("do_size_log",&value);
  switches->do_size_log=CONTENTS(value);

  get_named_variable("lumi_ee_2_n",&value);
  switches->lumi_ee_2_n=CONTENTS(value);

  get_named_variable("lumi_ee_2_min",&value);
  switches->lumi_ee_2_xmin=CONTENTS(value);

  get_named_variable("lumi_ee_2_max",&value);
  switches->lumi_ee_2_xmax=CONTENTS(value);

  get_named_variable("do_lumi_eg_2",&value);
  switches->do_lumi_eg_2=CONTENTS(value);

  get_named_variable("lumi_eg_2_n",&value);
  switches->lumi_eg_2_n=CONTENTS(value);

  get_named_variable("lumi_eg_2_min",&value);
  switches->lumi_eg_2_xmin=CONTENTS(value);

  get_named_variable("lumi_eg_2_max",&value);
  switches->lumi_eg_2_xmax=CONTENTS(value);

  get_named_variable("do_lumi_ge_2",&value);
  switches->do_lumi_ge_2=CONTENTS(value);

  get_named_variable("lumi_ge_2_n",&value);
  switches->lumi_ge_2_n=CONTENTS(value);

  get_named_variable("lumi_ge_2_min",&value);
  switches->lumi_ge_2_xmin=CONTENTS(value);

  get_named_variable("lumi_ge_2_max",&value);
  switches->lumi_ge_2_xmax=CONTENTS(value);

  get_named_variable("do_lumi_gg_2",&value);
  switches->do_lumi_gg_2=CONTENTS(value);

  get_named_variable("lumi_gg_2_n",&value);
  switches->lumi_gg_2_n=CONTENTS(value);

  get_named_variable("lumi_gg_2_min",&value);
  switches->lumi_gg_2_xmin=CONTENTS(value);

  get_named_variable("lumi_gg_2_max",&value);
  switches->lumi_gg_2_xmax=CONTENTS(value);

  get_named_variable("do_cross",&value);
  switches->do_cross=CONTENTS(value);

  get_named_variable("do_prod",&value);
  switches->do_prod=CONTENTS(value);

  get_named_variable("load_events",&value);
  switches->load_event=CONTENTS(value);

  get_named_variable("prod_e",&value);
  switches->prod_e=CONTENTS(value);

  get_named_variable("prod_scal",&value);
  switches->prod_scal=CONTENTS(value);

  get_named_variable("do_compt",&value);
  switches->do_compt=CONTENTS(value);

  get_named_variable("do_compt_phot",&value);
  switches->do_compt_phot=CONTENTS(value);

  get_named_variable("compt_x_min",&value);
  switches->compt_x_min=CONTENTS(value);

  get_named_variable("compt_scale",&value);
  switches->compt_scale=CONTENTS(value);

  get_named_variable("compt_emax",&value);
  switches->compt_emax=CONTENTS(value);

  get_named_variable("do_isr",&value);
  switches->do_isr=CONTENTS(value);

  get_named_variable("do_espread",&value);
  switches->do_espread=CONTENTS(value);

  get_named_variable("espread.1",&value);
  switches->espread1=CONTENTS(value);

  get_named_variable("which_espread.1",&value);
  switches->which_espread1=CONTENTS(value);

  get_named_variable("espread.2",&value);
  switches->espread2=CONTENTS(value);

  get_named_variable("which_espread.2",&value);
  switches->which_espread2=CONTENTS(value);

  get_named_variable("num_lumi",&value);
  switches->num_lumi=CONTENTS(value);

  get_named_variable("num_lumi_eg",&value);
  switches->num_lumi_eg=CONTENTS(value);

  get_named_variable("num_lumi_gg",&value);
  switches->num_lumi_gg=CONTENTS(value);

  get_named_variable("lumi_p",&value);
  switches->lumi_p=CONTENTS(value);

  get_named_variable("lumi_p_eg",&value);
  switches->lumi_p_eg=CONTENTS(value);

  get_named_variable("lumi_p_gg",&value);
  switches->lumi_p_gg=CONTENTS(value);

  get_named_variable("photon_ratio",&value);
  switches->photon_ratio=CONTENTS(value);

  get_named_variable("pair_ratio",&value);
  switches->pair_ratio=CONTENTS(value);

  get_named_variable("jet_ratio",&value);
  switches->jet_ratio=CONTENTS(value);

  get_named_variable("jet_ptmin",&value);
  switches->jet_pstar=CONTENTS(value);

  get_named_variable("store_jets",&value);
  switches->jet_store=CONTENTS(value);

  get_named_variable("jet_log",&value);
  switches->jet_select=CONTENTS(value);

  get_named_variable("beam_size",&value);
  switches->geom=CONTENTS(value);

  get_named_variable("beam_size_scale",&value);
  switches->r_scal=CONTENTS(value);

  get_named_variable("ext_field",&value);
  switches->ext_field=CONTENTS(value);

#ifndef EXT_FIELD
  if (switches->ext_field){
    fprintf(stderr,"Macro EXT_FIELD not defined\n");
    fprintf(stderr,"Cannot use flag ext_field\n");
    fprintf(stderr,"Please recompile defining Macro EXT_FIELD\n");
    exit(1);
  }
#endif

  get_named_variable("pair_q2",&value);
  switches->pair_q2=CONTENTS(value);

  get_named_variable("store_beam",&value);
  switches->store_beam=CONTENTS(value);

  get_named_variable("load_beam",&value);
  switches->load_beam=CONTENTS(value);

  get_named_variable("load_photons",&value);
  switches->load_photon=CONTENTS(value);

  get_named_variable("jet_pythia",&value);
  switches->jet_pythia=CONTENTS(value);

  get_named_variable("force_symmetric",&value);
  switches->force_symmetric=CONTENTS(value);

  get_named_variable("charge_symmetric",&value);
  switches->charge_symmetric=CONTENTS(value);

  get_named_variable("beam_vx_min",&value);
  switches->beam_vx_min=CONTENTS(value);

  get_named_variable("beam_vx_max",&value);
  switches->beam_vx_max=CONTENTS(value);

  get_named_variable("beam_vx_interval",&value);
  switches->beam_vx_interval=CONTENTS(value);

  get_named_variable("beam_vy_min",&value);
  switches->beam_vy_min=CONTENTS(value);

  get_named_variable("beam_vy_max",&value);
  switches->beam_vy_max=CONTENTS(value);

  get_named_variable("beam_vy_interval",&value);
  switches->beam_vy_interval=CONTENTS(value);

  get_named_variable("do_dump",&value);
  switches->do_dump=CONTENTS(value);

  get_named_variable("dump_step",&value);
  switches->dump_step=CONTENTS(value);

  get_named_variable("dump_particle",&value);
  switches->dump_particle=CONTENTS(value);

return 1;
}

#ifdef TWOBEAM
int file_read_parameters(DATEI *datei,char *name,BEAM_PARAMETERS *beam1,
                         BEAM_PARAMETERS *beam2,BEAM_PARAMETERS *beam1_2,
                         BEAM_PARAMETERS *beam2_2,GRID *grid,
                         SWITCHES *switches)
#else
int file_read_parameters(DATEI *datei,char *name,BEAM_PARAMETERS *beam1,
                         BEAM_PARAMETERS *beam2,GRID *grid,SWITCHES *switches)
#endif
{
  const int n=max_buffer2;
  char buffer[max_buffer2];
  float x;
  int m;
  do{	
      do{
          if (file_find_word(datei,"$PARAMETERS::")==0) return 0;
        }
      while(file_next_word(datei,name)==0);
   } while(!file_read_braces(datei,"{","}",buffer,n));
  def_param();
  input_values(buffer);
#ifdef TWOBEAM
  file_read_accelerator_i(beam1,beam2,beam1_2,beam2_2,grid,switches);
#else
  file_read_accelerator_i(beam1,beam2,grid,switches);
#endif
/*
  acc_test(beam1->ebeam,&(beam1->em_x),&(beam1->beta_x),&(beam1->sigma_x));
  acc_test(beam2->ebeam,&(beam2->em_x),&(beam2->beta_x),&(beam2->sigma_x));
  acc_test(beam1->ebeam,&(beam1->em_y),&(beam1->beta_y),&(beam1->sigma_y));
  acc_test(beam2->ebeam,&(beam2->em_y),&(beam2->beta_y),&(beam2->sigma_y));
*/
#ifdef TWOBEAM
  read_file_parameters_i(beam1,beam2,beam1_2,beam2_2,grid,switches);
#else
  read_file_parameters_i(beam1,beam2,grid,switches);
#endif
printf ("%g %g\n",beam1->beta_x,beam1->sigma_x);
printf ("%g %g\n",beam2->beta_x,beam2->sigma_x);
printf ("%g %g\n",beam1->beta_y,beam1->sigma_y);
printf ("%g %g\n",beam2->beta_y,beam2->sigma_y);
  return 1;
}

