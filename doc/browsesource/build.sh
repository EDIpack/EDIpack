#!/bin/bash


if [ -n "$1" ]; then
    RUNFORD=$1
else
    RUNFORD=false
fi


if $RUNFORD; then
cat <<EOF > ford_doc.md
---
project: EDIpack
preprocess:false
display:none
hide_undoc:true
print_creation_date:true
src_dir:../../src
output_dir:./ford_doc
extensions: f90
quiet:false
parallel:0
graph:true
graph_maxdepth:5
graph_maxnodes:20
---
This is my Fortran project!
EOF

rm -f src/*

#Run FORD (not checking actual presence)
ford ford_doc.md
#Sync module/*.html files in local src
rsync -avPhHO --del ford_doc/module/*.html src/
fi


rm -f ford_doc.md
rm -rf ford_doc

SRC=src
GRAPH=graphs
MOD=module
MOM=$(pwd)

rm -r $GRAPH
rm -r $MOD

mkdir -p $GRAPH
mkdir -p $MOD
rm $GRAPH/*
rm $MOD/*


#Generate the  .html images by stripping the svg part from the src files
cd $SRC
pwd 
for file in *.html;do
    echo $file
    sed -n "/<svg id/,/<\/svg>/p" $file > $MOM/$GRAPH/$file
done
cd $MOM


cd $GRAPH

#color tweaks:

#satisfy my ocd
for i in *.html; do
  sed -e "s/white/\#fcfcfc/g" $i > $i.new
  mv $i.new $i
done

#generic blocks in purple
for ifile in *.html; do
  #if there's link
  awk '/<title>/ {print; next_line=1; next} next_line {gsub(/#337ab7/, "#c061cb"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#337ab7/, "#c061cb"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done

#edipack in consistent blue. Use both regex
for ifile in *.html; do
  #if there's link
  awk '/title="ED/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/title="ED/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's link
  awk '/<title>ED/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>ED/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done

#edipack2ineq in consistent blue. Use both regex
for ifile in *.html; do
  #if there's link
  awk '/title="E2I/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/title="E2I/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's link
  awk '/<title>E2I/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#2980b9"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>E2I/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#2980b9"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done


#scifortran blocks in green 
for ifile in *.html; do
  #if there's link
  awk '/<title>SF_/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#26a269"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>SF_/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#26a269"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's link
  awk '/<title>scifor/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#26a269"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>scifor/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#26a269"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's link
  awk '/<title>SCIFOR/ {print; next_line=1; next} next_line {gsub(/#c061cb/, "#26a269"); next_line=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
  #if there's no link
  awk '/<title>SCIFOR/ {flag=2; print; next} flag == 2 {flag--; print; next} flag == 1 {gsub(/#c061cb/, "#26a269"); flag=0} 1' $ifile > $ifile.new
  mv $ifile.new $ifile
done

#make scifortran links work:
for ifile in *.html; do
  #SF modules
  awk '
  {
      print;
      if (match($0, /<g id="([^"]+)"/, arr)) {
          id_value = arr[1];
      }
      if (match($0, /<title>(SF_[^<]+)<\/title>/, arr)) {
          mytitle = arr[1];
          print "<g id=a_" id_value "><a xlink:href=\"https://scifortran.github.io/SciFortran/modules/" mytitle ".html\" xlink:title=\"" mytitle "\">";
      }
  }' $ifile > $ifile.new
  mv $ifile.new $ifile
  #scifor
  awk '
  {
      print;
      if (match($0, /<g id="([^"]+)"/, arr)) {
          id_value = arr[1];
      }
      if (match($0, /<title>scifor<\/title>/, arr)) {
          mytitle = arr[1];
          print "<g id=a_" id_value "><a xlink:href=\"https://scifortran.github.io/SciFortran/index.html\" xlink:title=\"" mytitle "\">";
      }
  }' $ifile > $ifile.new
  mv $ifile.new $ifile
  #SCIFOR
  awk '
  {
      print;
      if (match($0, /<g id="([^"]+)"/, arr)) {
          id_value = arr[1];
      }
      if (match($0, /<title>SCIFOR<\/title>/, arr)) {
          mytitle = arr[1];
          print "<g id=a_" id_value "><a xlink:href=\"https://scifortran.github.io/SciFortran/index.html\" xlink:title=\"" mytitle "\">";
      }
  }' $ifile > $ifile.new
  mv $ifile.new $ifile
done


#Generete a list of the actual .html files used in the src images, sorted and uniq-ed
grep  -ri "../module/" * | awk -F/ '/"/{ print $3 }' | awk '{print $1}' | sed -e "s/\"//g" |sort -u > $MOM/list2
cd $MOM
sed -e "s/html/rst/g" list2 > list
rm list2





#generate .rst

current_branch=$(git rev-parse --abbrev-ref HEAD)

for file in $(cat list); do
  name=$(echo $file | sed -e "s/\.rst//g")
  name_upper=$(echo $name | tr '[:lower:]' '[:upper:]' | sed -e "s/HXV/HxV/g")
  relativepath=$(find ../../ -type f -name "*${name_upper}.*90" | awk -Fsrc '{print $2}')
  if [ $name == "ed_version" ]; then
    relativepath=$(find ../../ -type f -name "*revision.in" | awk -Fsrc '{print $2}')  
  fi
  githubpath="https://github.com/EDIpack/EDIpack/tree/${current_branch}/src${relativepath}"
  echo $name_upper > module/$file
  echo "=====================================" >> module/$file
  echo " " >> module/$file
  if [ -f "graphs/${name}.html" ]; then
    echo "Found graph for $name"
    echo ".. raw:: html" >> module/$file
    echo "   :file:  ../graphs/${name}.html" >> module/$file
  echo " " >> module/$file
    echo "|" >> module/$file
  else
    echo "Not found graph for $name"
  fi
  echo " " >> module/$file
  echo "\`Open source file <${githubpath}>\`_ for :f:mod:\`${name}\` on GitHub" >> module/$file
  echo " " >> module/$file
done

# rm list
