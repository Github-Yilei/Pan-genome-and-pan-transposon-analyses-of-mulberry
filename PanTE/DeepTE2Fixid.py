import re
with open("opt_DeepTE.fasta", 'r') as lines:
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            new_line = re.sub(r"#.*__Class", "#Class", line)
            new_line_spl = new_line.split("#")
            # checking by hand: grep '>' | sed 's#.*\###' | sort | uniq
            ## class I
            if new_line_spl[1] == "ClassI":
                print(new_line_spl[0]  + "#" + "Unknown")
            elif new_line_spl[1] == "ClassI_LTR":
                print(new_line_spl[0]  + "#" + "LTR/unknown")
            elif new_line_spl[1] == "ClassI_LTR_Copia":
                print(new_line_spl[0]  + "#" + "LTR/Copia")
            elif new_line_spl[1] == "ClassI_LTR_Gypsy":
                print(new_line_spl[0]  + "#" + "LTR/Gypsy")
            elif new_line_spl[1] == "ClassI_nLTR":
                print(new_line_spl[0]  + "#" + "nonLTR/unknown")
            elif new_line_spl[1] == "ClassI_nLTR_DIRS":
                print(new_line_spl[0]  + "#" + "DIRS/unknown")
            elif new_line_spl[1] == "ClassI_nLTR_LINE_I":
                print(new_line_spl[0]  + "#" + "LINE/I")
            elif new_line_spl[1] == "ClassI_nLTR_LINE_L1":
                print(new_line_spl[0]  + "#" + "LINE/L1")
            elif new_line_spl[1] == "ClassI_nLTR_PLE":
                print(new_line_spl[0]  + "#" + "PLE/Penelope")
            elif new_line_spl[1] == "ClassI_nLTR_SINE_tRNA":
                print(new_line_spl[0]  + "#" + "SINE/tRNA")
            elif new_line_spl[1] == "LINE/L1__unknown":
                print(new_line_spl[0]  + "#" + "LINE/L1")
            elif new_line_spl[1] == "LTR/Caulimovirus__unknown":
                print(new_line_spl[0]  + "#" + "LTR/Caulimovirus")
            elif new_line_spl[1] == "LTR/Copia__unknown":
                print(new_line_spl[0]  + "#" + "LTR/Copia")
            elif new_line_spl[1] == "LTR/Gypsy__unknown":
                print(new_line_spl[0]  + "#" + "LTR/Gypsy")
            elif new_line_spl[1] == "LTR/unknown__unknown":
                print(new_line_spl[0]  + "#" + "LTR/unknown")
            ## class II
            elif new_line_spl[1] == "ClassII":
                print(new_line_spl[0]  + "#" + "Unknown")
            elif new_line_spl[1] == "ClassII_DNA_CACTA_MITE":
                print(new_line_spl[0]  + "#" + "TIR/CACTA")
            elif new_line_spl[1] == "ClassII_DNA_CACTA_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/CACTA")
            elif new_line_spl[1] == "ClassII_DNA_Harbinger_MITE":
                print(new_line_spl[0]  + "#" + "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "ClassII_DNA_Harbinger_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "ClassII_DNA_Harbinger_unknown":
                print(new_line_spl[0]  + "#" + "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "ClassII_DNA_hAT_MITE":
                print(new_line_spl[0]  + "#" + "TIR/hAT")
            elif new_line_spl[1] == "ClassII_DNA_hAT_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/hAT")
            elif new_line_spl[1] == "ClassII_DNA_hAT_unknown":
                print(new_line_spl[0]  + "#" + "TIR/hAT")
            elif new_line_spl[1] == "ClassII_DNA_Mutator_MITE":
                print(new_line_spl[0]  + "#" + "TIR/Mutator")
            elif new_line_spl[1] == "ClassII_DNA_Mutator_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/Mutator")
            elif new_line_spl[1] == "ClassII_DNA_Mutator_unknown":
                print(new_line_spl[0]  + "#" + "TIR/Mutator")
            elif new_line_spl[1] == "ClassII_DNA_P_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/P")
            elif new_line_spl[1] == "ClassII_DNA_TcMar_MITE":
                print(new_line_spl[0]  + "#" + "TIR/Tc1-mariner")
            elif new_line_spl[1] == "ClassII_DNA_TcMar_nMITE":
                print(new_line_spl[0]  + "#" + "TIR/Tc1-mariner")
            elif new_line_spl[1] == "ClassII_DNA_TcMar_unknown":
                print(new_line_spl[0]  + "#" + "TIR/Tc1-mariner")
            elif new_line_spl[1] == "ClassIII_Helitron":
                print(new_line_spl[0]  + "#" + "Helitron/helitron")
            elif new_line_spl[1] == "ClassII_MITE":
                print(new_line_spl[0]  + "#" + "MITE")
            elif new_line_spl[1] == "ClassII_nMITE":
                print(new_line_spl[0]  + "#" + "nonMITE")
            elif new_line_spl[1] == "DNA/CMC-EnSpm__unknown":
                print(new_line_spl[0]  + "#" + "CMC-EnSpm")
            elif new_line_spl[1] == "DNA/DTA__unknown":
                print(new_line_spl[0]  + "#" + "TIR/hAT")
            elif new_line_spl[1] == "DNA/DTC__unknown":
                print(new_line_spl[0]  + "#" + "TIR/CACTA")
            elif new_line_spl[1] == "DNA/DTH__unknown":
                print(new_line_spl[0]  + "#" + "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "DNA/DTM__unknown":
                print(new_line_spl[0]  + "#" + "TIR/Mutator")
            elif new_line_spl[1] == "DNA/DTT__unknown":
                print(new_line_spl[0]  + "#" + "TIR/Tc1-mariner")
            elif new_line_spl[1] == "DNA/Helitron__unknown":
                print(new_line_spl[0]  + "#" + "Helitron/helitron")
            elif new_line_spl[1] == "DNA/MULE-MuDR__unknown":
                print(new_line_spl[0]  + "#" + "MULE-MuDR")
            elif new_line_spl[1] == "DNA/PIF-Harbinger__unknown":
                print(new_line_spl[0]  + "#" + "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "MITE/DTA__unknown":
                print(new_line_spl[0]  + "#" +  "TIR/hAT")
            elif new_line_spl[1] == "MITE/DTC__unknown":
                print(new_line_spl[0]  + "#" +  "TIR/CACTA")
            elif new_line_spl[1] == "MITE/DTH__unknown":
                print(new_line_spl[0]  + "#" +  "TIR/PIF-Harbinger")
            elif new_line_spl[1] == "MITE/DTM__unknown":
                print(new_line_spl[0]  + "#" +  "TIR/Mutator")
            elif new_line_spl[1] == "RC/Helitron__unknown":
                print(new_line_spl[0]  + "#" +  "Helitron/helitron")
            ## Unknown
            elif new_line_spl[1] == "unknown__unknown":
                print(new_line_spl[0]  + "#" +  "Unknown")


        else:
            print(line)

#grep '>' panEDTA.TElib.Fixedid3ed.fa | sed 's#.*\###' | sort | uniq
