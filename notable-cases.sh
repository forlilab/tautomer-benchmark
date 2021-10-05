./tautomerizer.py -s "C1=CC=CN2C1=NC(=NC(=O)C)N2" -r pt07mod2.txt
./tautomerizer.py -s "N=C1c2n3c(-c4ccccc4C3)nc2NS(=O)(=O)N1" -r pt07mod2.txt
./tautomerizer.py -r pt07mod2.txt -s "NC1=c2c(nc3n2Cc2ccccc2-3)=NS(=O)(=O)N1"
N1N=CN2N=CN=C2-1
"COc1ccc(/N=N/c2c(C)c(C#N)c(=O)n(C)c2O)c([N+](=O)[O-])c1"
./tautomerizer.py -r pt07mod2.txt -s "OC1=CC(C)=NN1C2=CC=CC=C2"
./tautomerizer.py -r pt07mod2.txt -s "O=CC=NNC"
./tautomerizer.py -r pt07mod3.txt -s "CC1=NNC2N1N=CN=2" # example of PT_08_00

# PT_06 converts "CN1CCN2CC=C(O)N=C12" onto "CN1CCN2CCC(=O)N=C12" but not the reverse, by design!
