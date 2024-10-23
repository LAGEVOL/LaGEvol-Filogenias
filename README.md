# Processamento de Dados Genômicos e Inferências Filogenéticas
Repositório sobre o processamento de dados genômicos e inferências filogenéticas com dados gerados por sequenciamento de regiões-alvo com bibliotecas enriquecidas (*target-enrichment sequencing*). Todos os programas e comandos necessários devem ser executados no sistema operacional Linux (ou no Subsistema Linux para Windows), através da interface *Bash*.

Os comandos apresentados neste repositório foram baseados no workshop "Análises Filogenéticas Computacionais", ministrado por M. Kohler, M. Telhe e M. Romeiro-Brito, no Laboratório de Diversidade Genética e Evolução, UFSCar campus Sorocaba. 

## Índice

1. [Trimagem inicial dos dados brutos](#Trimagem-inicial-dos-dados-brutos)

2. [Montagem dos dados trimados em locos](#Montagem-dos-dados-trimados-em-locos)

3. [Identificação e remoção de locos possivelmente parálogos](#Identificação-e-remoção-de-locos-possivelmente-parálogos)

4. [Alinhamento de sequências](#Alinhamento-de-sequências)

5. [Polimento das sequências alinhadas](#Polimento-das-sequências-alinhadas)

6. [Inferência filogenética utilizando o método de máxima verossimilhança](#Inferência-filogenética-utilizando-o-método-de-máxima-verossimilhança)

   6.1 [Remoção de regiões hipervariáveis](#Remoção-de-regiões-hipervariáveis)

7. [Inferência filogenética utilizando o método de coalescência](#Inferência-filogenética-utilizando-o-método-de-coalescência)

## 1. Trimagem inicial dos dados brutos
Definindo as variáveis referentes às duas fitas do DNA sequenciado:

```
R1=(*R1.fastq)
R2=(*R2.fastq)
```

Para esta trimagem inicial dos dados brutos, utilizaremos o *fastp*. Este programa remove sequências de baixa qualidade, sequências muito curtas e adaptadores dos nossos dados. Estes parâmetros podem ser modificados com comandos específicos, mas neste caso utilizaremos os valores padrões fornecidos pelo programa. Para mais informações, acesse: https://github.com/OpenGene/fastp.git.

Utilizando o *fastp* para uma única amostra:

```
fastp -i /caminho/para/a/pasta -I /caminho/para/a/pasta -o /caminho/para/a/pasta/_trimmed.fastq -O /caminho/para/a/pasta/_trimmed.fastq
```

Utilizando o *fastp* para várias amostras: 

```
for ((i=0;i<=${#R1[@]};i++)); do fastp  -i "${R1[i]}" -I "${R2[i]}"  -o "${R1[i]%.fastq}_trimmed.fastq" -O "${R2[i]%.fastq}_trimmed.fastq" -j "${R1[i]}.fastp.json" -h "${R1[i]}.fastp.html" -q 20 --dont_overwrite --failed_out "failed.${R1[i]}"; done
```

Os arquivos com os dados trimados terão o sufixo `_trimmed.fastq`. Para melhor organização, é recomendável que estes arquivos sejam movidos para uma nova pasta. 

## 2. Montagem dos dados trimados em locos
Para montar os dados trimados em seus respectivos locos (processo chamado de *assembling*), utilizaremos o programa *HybPiper* (https://github.com/mossmatters/HybPiper.git). Precisaremos de um arquivo de referência para a montagem dos locos, indicado nos comandos pela *flag* `-t_dna`. 

Utilizando o *HybPiper* em uma única amostra:

```
hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix nome_da_amostra --bwa --no_padding_supercontigs --thresh 75
```

Utilizando o *HybPiper* em várias amostras:

Primeiramente, crie um arquivo listando as suas amostras (ex: namelist.txt).
```
while read name; do hybpiper assemble -r amostra_trimada_R1.fastq amostra_trimada_R2.fastq -t_dna arquivo_de_referência.fasta --prefix $name --bwa --no_padding_supercontigs --thresh 75; done < namelist.txt
```

O *HybPiper* também possui comandos para gerar estatísticas sobre os locos montados.

Gerando estatísticas no formato de tabela (seq_lengths.tsv):
```
hybpiper stats -t_dna arquivo_de_referência.fasta gene namelist.txt
```
Gerando estatísticas no formato de mapa de calor:
```
hybpiper recovery_heatmap seq_lengths.tsv
```

Os *reads* recuperados para cada região alvo devem ser sintetizados em um único arquivo multi-fasta, com o comando:

```
hybpiper retrieve_sequences dna -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_regioes_alvo
```

Utilizando o *HybPiper*, também é possível recuperar regiões adjacentes às regiões alvo do painel de captura utilizado no sequenciamento. A combinação de regiões alvo do painel + sequências adjacentes é chamada de *supercontig*, e pode ser de grande valor a depender dos dados e análises trabalhadas.
```
hybpiper retrieve_sequences supercontig -t_dna arquivo_de_referência.fasta --sample_names namelist.txt --fasta_dir nome_da_pasta_supercontigs
```

Analisando as estatísticas geradas pelo *HybPiper*, você pode encontrar amostras que não atingiram os parâmetros esperados. Neste caso, você pode removê-las da sua base de dados:

Crie um arquivo listando as amostras que deseja remover (ex: samples_to_delete.txt).
```
for i in *.fasta; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
for i in *.FNA; do pxrms -s $i -f /caminho/para/a/pasta/samples_to_delete.txt -o "${i%.FNA}_Reduced.fasta"; done
```

Por fim, é importante que você remova avisos que o programa pode ter adicionado ao nome das suas amostras. Para isso, utilize os seguintes comandos para os dois formatos de arquivo gerados (fasta e FNA):
```
sed -i 's/\ single_hit//g' *.fasta
sed -i -E 's/\ multi_hit_stitched_contig_comprising_[0-9]+_hits//g' *.fasta
sed -i 's/\.0.*/.0/g' *.fasta
sed -i 's/\.main.*/.main/g' *.fasta
sed -i 's/ .*//' *.fasta
sed -i 's/\./@/g' *.fasta
```

## 3. Identificação e remoção de locos possivelmente parálogos
A identificação de locos possivelmente parálogos pode ser realizada de diversas maneiras. No presente repositório, adotaremos a análise a olho das árvores de gene, buscando por indícios de paralogia em cada um dos locos. Para isto, precisamos utilizar o *HybPiper* para recuperar sequências que o programa identificou como parálogas durante o *assembling*:
```
hybpiper paralog_retriever namelist.txt -t_dna arquivo_de_referência.fasta
```
Este comando irá gerar estatísticas sobre os parálogos detectados pelo programa na sua base de dados, além de uma pasta contendo todas as sequências identificadas como parálogas, nomeada como `paralogs_all`. As sequências desta pasta podem ser alinhadas (com o *MAFFT*) e usadas para gerarmos árvores de cada um dos locos (com o *IQTree*):
```
nohup sh -c 'for i in *.fasta; do mafft --reorder --auto "$i" > "caminho/para/a/pasta/aligned_$i"; done'  &
iqtree -s aligned_nome_da_amostra_paralogs_all.fasta
```
As árvores geradas devem ser analisadas manualmente em busca de indícios de paralogia, como múltiplas cópias de uma amostra com posicionamento distante na mesma árvore. Para mais detalhes sobre a identificação a olho de parálogos, acesse: https://doi.org/10.1093/sysbio/syad076. 

![image](https://github.com/user-attachments/assets/d74a76e8-4f6c-4b0c-92f7-2219696fd03a)

Após a identificação dos locos possivelmente parálogos, você pode removê-los da sua base de dados, movendo-os para uma pasta distinta:

Crie uma lista com os locos que você identificou como possivelmente parálogos (ex: paralogs_list.txt).
```
mkdir locos_paralogos
while read line; do mv $line ./locos_paralogos; done < paralogs_list.txt
```

## 4. Alinhamento de sequências
Para alinhar as sequências de um determinado loco obtidas para cada uma das amostras, utilizaremos o programa *MAFFT* (https://github.com/GSLBiotech/mafft). É recomendável que você crie novas pastas para armazenar os alinhamentos. No presente repositório, utilizaremos os parâmetros padrões fornecidos pelo programa:

```
nohup sh -c 'for i in *.fasta; do mafft --reorder --auto "$i" > "caminho/para/a/pasta/aligned_$i"; done'  &
```
Os alinhamentos gerados terão o prefixo `aligned_`.

## 5. Polimento das sequências alinhadas
Uma segunda trimagem dos dados é necessária para remover posições dos alinhamentos com uma alta proporção de *gaps* e para isso, utilizaremos o *trimal* (https://github.com/inab/trimal.git). A *flag* que determina o limite de *gaps* a ser tolerado é a `-gt`. Neste caso, `-gt 0.7` indica que serão removidas as "colunas" do alinhamento que possuem uma fração de *gaps* maior ou igual a 30%:

```
mkdir alinhamento_trimado
for i in *.fasta; do trimal -in $i -out ./alinhamento_trimado/"$i"_trimmed.fasta -gt 0.7; done;
```

O AMAS (https://github.com/marekborowiec/AMAS.git) pode ser utilizado para gerar estatísticas dos seus dados trimados, na forma de uma planilha. Estas estatísticas podem ser úteis caso você deseje explorar o efeito dos dados faltantes (*gaps*) sobre a sua base de dados.

```
python3 AMAS.py summary -f fasta -d dna -i *.fasta -o SummaryStats.csv
```

## 6. Inferência filogenética utilizando o método de máxima verossimilhança
Inicialmente, o programa *IQtree* (https://github.com/iqtree/iqtree2.git) será utilizado para gerar uma árvore para cada loco:

```
nohup sh -c 'for i in *.fasta; do iqtree2 -nt 4 -s "$i" -st DNA -m MFP -B 10000; done 2>iqtree.err' &
```

No método de máxima verossimilhança, as árvores geradas para cada loco precisam ser concatenadas, formando uma supermatriz:
```
pxcat -s *fasta -p partitions.txt -o minha_supermatriz.fasta
```
A *flag* `-p` gera um arquivo de partições, que pode ser útil caso você precise "desconcatenar" a supermatriz posteriormente.

Uma árvore de espécies pode ser gerada a partir da supermatriz:
```
iqtree -nt 4 -s minha_supermatriz.fasta -st DNA -m MFP -B 10000
```
### 6.1 Remoção de regiões hipervariáveis
É recomendável que uma última trimagem dos dados seja aplicada, utilizando o *spruceup* (https://github.com/marekborowiec/spruceup.git). Este programa remove regiões hipervariáveis (*outliers*) da nossa base de dados. Para utilizar o programa precisaremos da supermatriz, de uma árvore de espécies (opcional) e de um arquivo de texto com os parâmetros a serem adotados. Um modelo deste último arquivo está anexado ao repositório (`configuration_spruceup.conf`).

```
python -m spruceup configuration-spruceup.conf
```

Este comando irá gerar uma supermatriz trimada (`0.valor_do_cut_off_minha_supermatriz.fasta`) que pode ser utilizada diretamente para gerar uma nova árvore de espécies, pelo método de máxima verossimilhança. A supermatriz trimada também pode ser "desconcatenada" utilizando o arquivo de partições (`partitions.txt`) gerado anteriormente:
```
python3 AMAS.py split -f fasta -d dna -i supermatriz_trimada.fasta -l partitions.txt -u fasta -j
```
Esta etapa é importante para que você possa inferir árvores de espécie pelo método de coalescência. 
O AMAS também pode ser utilizado para gerar estatísticas nesta etapa:
```
python3 AMAS.py summary -f fasta -d dna -i *.fasta -o SummaryStats.csv
```

## 7. Inferência filogenética utilizando o método de coalescência
Para gerar uma árvore de espécies coalescente, precisamos reunir todas as árvores de gene em um único arquivo (all_gene.trees):
```
cat *.treefile > all_gene.trees
```
É recomendável que este arquivo contendo todas as árvores reunidas seja polido, para reduzir a influência de nós pouco informativos na nossa árvore de espécies. Para isso, usaremos um comando que colapsa nós que não atingem um valor de suporte mínimo:
```
nw_ed all_gene.trees 'i & b<=50' o >  all_gene_BS50.tree
```
Neste comando foi estabelecido um suporte mínimo de 50 para os nós, mas o limite ideal pode variar de acordo com os nossos dados. Para mais informações, acesse: https://doi.org/10.1016/j.ympev.2021.107092.

Com as nossas árvores de gene unidas e polidas, utilizzaremos o ASTRAL III (https://github.com/smirarab/ASTRAL.git) para gerar nossa árvore de espécies coalescente:
```
java astral.5.7.8.jar -i all_gene_BS50.tree -o sptree_astral_BS50.tree 2> sptree_astral.log
```










