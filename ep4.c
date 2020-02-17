/* ************************************************************************** */
/* Nome: Giovana Oshiro da Silva                                              */
/* Numero USP: 8022103                                                        */
/*                                                                            */
/* Nome: Lucas Freitas Bastos                                                 */
/* Numero USP: 9783118                                                        */
/*                                                                            */
/* Exercicio-programa 4                                                       */
/* ************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#define MAX 2000303

struct reg{
  int score;
  struct reg *diagonal, *lado, *cima;
};

typedef struct reg celula; /* estrutura da matriz */

/* A funcao "similaridade" completa a matriz de pontuacao a partir do algoritmo
Needleman-Wunsch. Cada casa da matriz eh uma estrutura do tipo celula, que arma-
zena as informacoes de score e para qual casa (de cima, de baixo ou diagonal) uma
casa especifica da matriz aponta (a partir do algoritmo). */

void similaridade (celula **pontuacao, int pesos[20][20], int *seq1, int *seq2, int tam_seq1, int tam_seq2, int gap){
  int i, j=0, k, aux [3], max = 0;

  pontuacao[0][0].score = 0;

  for (i = 1; i < tam_seq1; i++){
    pontuacao[i][0].score = i*gap;
    pontuacao[i][0].cima = &pontuacao[i-1][0];
  }

  for (j = 1; j < tam_seq2; j++){
    pontuacao[0][j].score = j*gap;
    pontuacao[0][j].lado = &pontuacao[0][j-1];
  }

  for (i = 1; i < tam_seq1; i++){
    for (j = 1; j < tam_seq2; j++){
      aux [0] = pontuacao[i-1][j].score + gap;
      aux [1] = pontuacao[i-1][j-1].score + pesos [seq1[i]][seq2[j]];
      aux [2] = pontuacao[i][j-1].score + gap;
      max = aux [0];
      for (k = 1; k < 3; k++) if (aux [k] > max) max = aux [k];
      pontuacao[i][j].score = max;
      if (aux[0] == max) pontuacao[i][j].cima = &pontuacao[i-1][j];
      if (aux[1] == max) pontuacao[i][j].diagonal = &pontuacao[i-1][j-1];
      if (aux[2] == max) pontuacao[i][j].lado = &pontuacao[i][j-1];
    }
  }
}

/* A funcao "caminhos" retorna o numero total de melhores alinhamentos encontrados.
Ela percorre a matriz a partir da ultima posicao e verifica se os ponteiros "cima",
"lado" e "diagonal" sao nulos ou nao. Caso o primeiro ponteiro verficado (no nosso
codigo, esse ponteiro eh "cima") nao seja nulo, a funcao eh chamada recursivamente
com a casa para qual esse ponteiro esta apontando e assim por diante. Desse modo,
a funcao percorre todos os caminhos que levam da ultima casa da matriz ate a pri-
meira. Conforme o ponteiro "anda", as correspondencias (match, mismatch ou gap)
sao analisadas e armazenadas em 4 vetores, para que depois sejam impressos. Quando
todos os ponteiros forem nulos, significa que um caminho foi encontrado e, portanto,
ele eh impresso. */

int caminhos (int i, int j, int k, celula *aux, char seq1 [], char seq2 [], char mismatch [], char print_seq1 [], char match [], char print_seq2 [], FILE * saida, int alinhamentos){

  if (aux->cima == NULL && aux->diagonal == NULL && aux->lado == NULL){ /* caminho eh impresso no formato FASTA */
    int t = 0, p, cont = k-1;
    fprintf (saida, "%dº Alinhamento:\n", alinhamentos);
    while (t >= 0){
      t = cont;
      for (p = 0; t >= 0 && p < 80; p++, t--){
        fprintf (saida, "%c", mismatch [t]);
      }
      t = cont;
      fprintf (saida, "\n");
      for (p = 0; t >= 0 && p < 80; p++, t--) fprintf (saida, "%c", print_seq1 [t]);
      t = cont;
      fprintf (saida, "\n");
      for (p = 0; t >= 0 && p < 80; p++, t--) fprintf (saida, "%c", match [t]);
      t = cont;
      fprintf (saida, "\n");
      for (p = 0; t >= 0 && p < 80; p++, t--) fprintf (saida, "%c", print_seq2 [t]);
      cont = t;
      fprintf (saida, "\n");
      fprintf (saida, "\n");
    }
    alinhamentos++;
    return alinhamentos;
  }

  if (aux->cima != NULL){ /* se aux->cima nao eh nulo, colocamos um gap na sequencia 2 */
    mismatch [k] = ' ';
    print_seq1 [k] = seq1 [i];
    match [k] = ' ';
    print_seq2 [k] = '-';
    alinhamentos = caminhos (i-1, j, k+1, aux->cima, seq1, seq2, mismatch, print_seq1, match, print_seq2, saida, alinhamentos);
  }

  if (aux->diagonal != NULL){ /* se aux->diagonal nao eh nulo, verificamos se temos
                                um match ou um mismatch */
    if (seq1 [i] == seq2 [j]){
      mismatch [k] = ' ';
      print_seq1 [k] = seq1 [i];
      match [k] = '|';
      print_seq2 [k] = seq2 [j];
    }
    else{
      mismatch [k] = '*';
      print_seq1 [k] = seq1 [i];
      match [k] = ' ';
      print_seq2 [k] = seq2 [j];
    }
    alinhamentos = caminhos (i-1, j-1, k+1, aux->diagonal, seq1, seq2, mismatch, print_seq1, match, print_seq2, saida, alinhamentos);
  }

  if (aux->lado != NULL){ /* se aux->lado nao eh nulo, inserimos um gap na sequencia 1 */
    mismatch [k] = ' ';
    print_seq1 [k] = '-';
    match [k] = ' ';
    print_seq2 [k] = seq2 [j];
    alinhamentos = caminhos (i, j-1, k+1, aux->lado, seq1, seq2, mismatch, print_seq1, match, print_seq2, saida, alinhamentos);
  }

  return alinhamentos;
}

int main (){

  FILE * entrada;
  FILE * saida;
  celula **pontuacao, *aux;
  char caractere, arquivo_entrada [50], arquivo_saida [50], *seq1, *seq2;
  char mismatch [MAX], print_seq1 [MAX], match [MAX], print_seq2 [MAX];
  int gap, *sequencia1, *sequencia2, leitura = 0, flag = 0, tam_seq1 = 1, tam_seq2 = 1;
  int i, j, k, alinhamentos = 1;
  int aminoacidos [27] = {-1, 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, -1, 18, -1};
      /* Completamos o vetor da seguinte forma: cada casa representa uma letra do
        alfabeto (em ordem crescente, isto eh, de A=1 em diante); se a casa-letra
        eh simbolo de um aminoacido, ela recebe o numero (a partir do 0) da posicao
        do aminoacido referente a tabela de pesos (que tambem eh a ordem alfabetica).
        Caso contrario, a casa recebe -1 (nao ha aminoacidos com aquele simbolo).
        Por exemplo, a casa 1 representa a letra A; A eh simbolo do aminoacido
        alanina, que eh o primeiro aminoacido em ordem alfabetica. Portanto, a
        casa 1 do vetor "aminoacidos" recebe 0. */
  int pesos [20][20] = /* matriz BLOSUM62 */
  {{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
   {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
   {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
   {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
   { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
   {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
   {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
   { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
   {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
   {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
   {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
   {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
   {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
   {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
   {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
   { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
   { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
   {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
   {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
   { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};

  printf ("Entre com o nome do arquivo que contem as sequencias no formato FASTA:\n");
  scanf ("%s", arquivo_entrada);

  printf ("Entre com o nome do arquivo de saida:\n");
  scanf ("%s", arquivo_saida);

  printf ("Entre com a penalidade para o gap:\n");
  scanf ("%d", &gap);

  entrada = fopen (arquivo_entrada, "r");
  saida = fopen (arquivo_saida, "w");

  /* Serao criados dois vetores para a sequencia 1 e dois vetores para a sequencia
  2. Em um deles a sequencia sera armazenada com seus codigos em char e, no outro,
  teremos a sua correspondencia em int. Esses vetores sao alocados dinamicamente. */

  sequencia1 = malloc (1*sizeof(int));
  sequencia2 = malloc (1*sizeof(int));
  seq1 = malloc (1*sizeof(char));
  seq2 = malloc (1*sizeof(char));

  sequencia1 [0] = -1;
  sequencia2 [0] = -1;

  /* Leitura do arquivo e armazenamento das sequencias nos vetores. */

  fscanf (entrada, "%c", &caractere);

  while (caractere != '\n') fscanf (entrada, "%c", &caractere);

  fscanf (entrada, "%c", &caractere);

  while (flag != 2){ /* armazenamento da sequencia 1 em vetores */
    if (caractere == '\n'){
      flag++;
      fscanf (entrada, "%c", &caractere);
    }
    else{
      flag = 0;
      sequencia1 = realloc (sequencia1, (tam_seq1+1)*sizeof(int));
      sequencia1 [tam_seq1] = aminoacidos [caractere-64];
      seq1 = realloc (seq1, (tam_seq1+1)*sizeof(char));
      seq1 [tam_seq1] = caractere;
      tam_seq1++;
      fscanf (entrada, "%c", &caractere);
    }
  }

  while (caractere != '\n') fscanf (entrada, "%c", &caractere);

  fscanf (entrada, "%c", &caractere);

  while (leitura != EOF){ /* armazemento sequencia 2 */
    if (caractere == '\n'){
      leitura = fscanf (entrada, "%c", &caractere);
    }
    else{
      sequencia2 = realloc (sequencia2, (tam_seq2+1)*sizeof(int));
      sequencia2 [tam_seq2] = aminoacidos [caractere-64];
      seq2 = realloc (seq2, (tam_seq1+1)*sizeof(char));
      seq2 [tam_seq2] = caractere;
      tam_seq2++;
      leitura = fscanf (entrada, "%c", &caractere);
    }
  }

  pontuacao = malloc (tam_seq1*sizeof(celula)); /* matriz alocada dinamicamente */

  for (i = 0; i < tam_seq1; i++){ /* matriz tam_seq1 x tam_seq2 criada */
    pontuacao [i] = malloc (tam_seq2*sizeof(celula));
    for (j = 0; j < tam_seq2; j++){
      pontuacao[i][j].score = 0;
      pontuacao[i][j].lado = NULL;
      pontuacao[i][j].cima = NULL;
      pontuacao[i][j].diagonal = NULL;
    }
  }

  /* Matriz eh preenchida a partir do algoritmo Needleman-Wunsch. */
  similaridade (pontuacao, pesos, sequencia1, sequencia2, tam_seq1, tam_seq2, gap);

  for (i = 0; i < MAX; i++){ /* "zerando" vetores de char */
    mismatch [i] = '\0';
    print_seq1 [i] = '\0';
    match [i] = '\0';
    print_seq2 [i] = '\0';
  }

  aux = &pontuacao[tam_seq1-1][tam_seq2-1];

  i = tam_seq1-1;
  j = tam_seq2-1;
  k = 0;

  fprintf (saida, "SCORE: %d\n", pontuacao[tam_seq1-1][tam_seq2-1].score);
  fprintf (saida, "\n");

  /* Caminhos sao percorridos e o numero total eh retornado. Dentro da funcao, os
  caminhos sao impressos. */
  alinhamentos = caminhos (i, j, k, aux, seq1, seq2, mismatch, print_seq1, match, print_seq2, saida, alinhamentos);

  fprintf (saida, "Foi encontrado um total de %d melhor(es) alinhamento(s).\n", alinhamentos-1);

  printf ("Verifique o arquivo de saida para visualizar os melhores alinhamentos.\n");

  fclose (entrada);
  fclose (saida);

  return 0;

}
