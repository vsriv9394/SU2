
cp m6.meshb tmp.meshb

amg -in tmp -geotol  0.0001  -hgrad 1.2 -hmax 4 -prepro  -surf 1 

#mv tmp.o.back.meshb CRM_back.meshb
#mv tmp.o.back.solb CRM_back.solb
#
#amg -in tmp -hmsh -cfacaniso 2.0 -back CRM_back -surf 3 -hmax 70 -hgrad 1.5 -nordg