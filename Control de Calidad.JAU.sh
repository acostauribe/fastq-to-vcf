#!/bin/bash

##SCRIPT DE CONTROL DE CALIDAD
#Natalia Acosta Baena y Juliana Acosta Uribe
#Marzo 2023

#===========
ssh acostabaena@kosik.cnsi.ucsb.edu
Federico.2023
ssh node10
#=============
#variables

wd='/home/acostabaena/QC'
fastqc='/home/acostabaena/bin/FastQC/fastqc'
trimmomatic='/home/acostabaena/bin/Trimmomatic-0.39/trimmomatic-0.39.jar'

# PREPARAR DATOS
## Mover archivos al directorio de Input
mkdir $wd/fastqc_input
mv *.fastq.gz $wd/fastqc_input

## Generar una lista de los Prefijos (todo lo que dice antes de fastq.gz). Cada prefijo representa una muestra
cd $wd/fastqc_input

ls > fastq_files.txt #generar una lista de todos los archivos. 
#Cada muestra esta representada por dos archivos  {muestra}_1.fastq.gz y {muestra}_2.fastq.gz
sed -i "s/_1.fastq.gz//g" fastq_files.txt #borrar el sufijo "_1.fastq.gz"
sed -i "s/_2.fastq.gz//g" fastq_files.txt #borrar el sufijo "_2.fastq.gz"
uniq fastq_files > tmp  #borrar duplicados
mv tmp fastq_files.txt #renombrar archivo

## Despues de este paso quedamos con una lista de todos los prefijos de los fastq para iterar

# PASO 1: Revisar calidad de los fastq con fastQC

## 1. Visualizacion de la calidad de los datos crudos con FastQC
## asumiendo que el archivo esta en $wd/fastqc_input_raw
cd $wd/fastqc_input_raw
mkdir $wd/fastqc_output_raw

## Utilizar la funcion 'while read' para leer todos los prefijos en fastq_files.txt [batch]
while read line; 
do 
$fastqc ${line}_1.fastq.gz ${line}_2.fastq.gz -o $wd/fastqc_output_raw ;
done < fastq_files.txt 

## Analizar individualmente 
#$fastqc 378F4AT_1.fastq.gz 378F4AT_2.fastq.gz -o $wd/fastqc_output_raw

## 'fastqc' genera dos archivos por muestra (line): {muestra}_1_fastqc.zip y {muestra}_1_fastqc.html
## Opcional: organizar output de fastqc
cd $wd/fastqc_output_raw 
mkdir zip
mv *.zip $wd/fastqc_output_raw/zip
mkdir html
mv *.html $wd/fastqc_output_raw/html


# PASO 2: Edicion de secuencias de mala calidad con TRIMMOMATIC

##Crear carpeta para Trimmomatic:
mkdir Trimmomatic
## Crear otra dentro de Trimmomatic para los output
mkdir trimmomatic_output 
cd $wd/Trimmomatic


## Batch  (procesar todos las muestras escritas en fastq_files.txt)
trimmomatic=/home/acostabaena/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
wd=/home/acostabaena/QC
while read line; 
do 
java -jar $trimmomatic PE -threads 20 -phred33 -trimlog ${line}_trimmomatic.log \
${wd}/fastqc_input_raw/${line}_1.fastq.gz ${wd}/fastqc_input_raw/${line}_2.fastq.gz \
${line}_1P.fastq.gz ${line}_1U.fastq.gz \
${line}_2P.fastq.gz ${line}_2U.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
done < $wd/fastqc_input_raw/fastq_files.txt 

##trimmomatic2
trimmomatic=/home/acostabaena/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
wd=/home/acostabaena/QC
while read line; 
do 
java -jar $trimmomatic PE -threads 20 -phred33 -trimlog ${line}_trimmomatic.log \
${wd}/fastqc_input/${line}_1.fastq.gz ${wd}/fastqc_input/${line}_2.fastq.gz \
${line}_1P.fastq.gz ${line}_1U.fastq.gz \
${line}_2P.fastq.gz ${line}_2U.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ;
done < $wd/fastqc_input/fastq_files.txt

## El input son las dos secuencias de una muestra.
## El output son 4 archivos, por cada lectura sale un (Prefijo_1)U.fastq.gz (Prefijo_2)U.fastq.gz(unpaired) y un _P.fastq.gz (paired)

## comando para pasar de una carpeta a otra:
mv *.gz /home/acostabaena/QC/Trimmomatic/trimmomatic_output
mv *.log /home/acostabaena/QC/Trimmomatic/trimmomatic_output


## Si se quiere procesar un numero pequeño de muestras, se puede hace con una 'for loop'. 
## Se puede poner una sola muestra y hace una analisis "individual"


# PASO 3: Revisar NUEVAMENTE LA CALIDAD de los fastq editados con trimmomatic usando fastQC 

## El input son las dos secuencias P (paired) de una muestra: _1P.fastq.gz _2P.fastq.gz 
## El output son los dos archivos de fastqc por muestra: {muestra}_1P_fastqc.zip y {muestra}_1P_fastqc.html

mkdir fastqc_output_trimmed

cd $wd/fastqc_output_trimmed

## Utilizar la funcion 'while read' para leer todos los prefijos en trimm_files.txt [batch]

wd=/home/acostabaena/QC
cd $wd/Trimmomatic/trimmomatic_output
while read line; 
do 
$fastqc ${line}_1P.fastq.gz ${line}_2P.fastq.gz -o $wd/fastqc_output_trimmed ;
done < $wd/fastqc_input_raw/fastq_files.txt 


wd=/home/acostabaena/QC
cd $wd/trimmomatic2
while read line; 
do 
$fastqc ${line}_1P.fastq.gz ${line}_2P.fastq.gz -o $wd/fastqc_output_trimmed2 ;
done < $wd/fastqc_input_raw/fastq_files.txt 

mkdir zip
mv *.zip $wd/fastqc_output_raw/zip
ls $wd/fastqc_output_raw/zip > fastqc_output_zip_files.txt
sed -i "s/.zip//g" fastqc_output_zip_files.txt #borrar el sufijo ".zip"

### descomprimir los zip y juntar los summaries en un solo archivo

while read line; 
do 
#unzip ${line}.zip
cp ./${line}/summary.txt ${line}_summary.txt
awk '{print $1}' ${line}_summary.txt > ${line}_summary_column1.txt
done < fastqc_output_zip_files.txt

paste *_summary_column1.txt > all_column1.txt
tr '\n' '\t' < fastqc_output_zip_files.txt > header.txt
cat header.txt all_column1.txt > all_files_summary_column1.txt 

## Analizar individualmente 
##$fastqc 1_1P.fastq.gz 1_2P.fastq.gz -o $wd/trimmomatic_output_trimmed