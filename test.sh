path="/mnt/d/desktop/"
in="titanic.csv"
out="out.csv"
filt="0,3,Mr."

#cat $path$in | parallel --pipe awk \'{print \$1}\' | awk '{print $1}'

cat $path$in | parallel --pipe -j 2 awk -va=$filt \'{if\(\$1 == a  \&\& \$2 == \"Owen\" \) print \$0}\' | awk -va=$filt '{if($1 == a && $2 == "Owen") print $0}' 
