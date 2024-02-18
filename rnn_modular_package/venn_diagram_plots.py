"""
module to plot venn diagrams.

Input file:
    output file from rank_and_comparison.py script in csv format

Limitations:
    Work only upto three trials
    works only with python2.7

"""
import os
import sys
import getopt
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
from collections import OrderedDict


def read_and_plot_venn_diagram(cname='', output_filename='output.png'):
    # read csv file
    df = pd.read_csv(cname, delimiter='\t')
    # print (df.columns)
    df_dict = OrderedDict()
    for a in df.columns:
        if 'rank' in a:
            print(a)
            df_dict.update({a: list(df[a])})
    if float(len(df.keys()) - 2) / 2 == 3:
        # calculate intersections
        
        ABc = 0
        AbC = 0
        aBC = 0
        ABC = 0
        
        for i in range(0, len(df_dict[df_dict.keys()[0]])):
            if df_dict[df_dict.keys()[0]][i] == df_dict[df_dict.keys()[1]][i] and df_dict[df_dict.keys()[0]][i] != \
                    df_dict[df_dict.keys()[2]][i] and df_dict[df_dict.keys()[1]][i] != df_dict[df_dict.keys()[2]][i]:
                ABc += 1
            
            if df_dict[df_dict.keys()[0]][i] != df_dict[df_dict.keys()[1]][i] and df_dict[df_dict.keys()[0]][i] == \
                    df_dict[df_dict.keys()[2]][i] and df_dict[df_dict.keys()[1]][i] != df_dict[df_dict.keys()[2]][i]:
                AbC += 1
            
            if df_dict[df_dict.keys()[0]][i] != df_dict[df_dict.keys()[1]][i] and df_dict[df_dict.keys()[0]][i] != \
                    df_dict[df_dict.keys()[2]][i] and df_dict[df_dict.keys()[1]][i] == df_dict[df_dict.keys()[2]][i]:
                aBC += 1
            
            if df_dict[df_dict.keys()[0]][i] == df_dict[df_dict.keys()[1]][i] and df_dict[df_dict.keys()[0]][i] == \
                    df_dict[df_dict.keys()[2]][i] and df_dict[df_dict.keys()[1]][i] == df_dict[df_dict.keys()[2]][i]:
                ABC += 1
        
        Abc = len(df_dict[df_dict.keys()[0]]) - ABc - AbC - ABC
        aBc = len(df_dict[df_dict.keys()[1]]) - ABc - aBC - ABC
        abC = len(df_dict[df_dict.keys()[2]]) - aBC - ABc - ABC
        
        #
        # feed parameters to create venn diagram
        #
        
        s = (
            Abc,
            aBc,
            ABc,
            abC,
            AbC,
            aBC,
            ABC,
        )
        
        v = venn3(subsets=s, set_labels=('A', 'B', 'C'))
        
        # Subset labels
        v.get_label_by_id('100').set_text('Abc ' + str(Abc))
        v.get_label_by_id('010').set_text('aBc ' + str(aBc))
        v.get_label_by_id('110').set_text('ABc ' + str(ABc))
        v.get_label_by_id('001').set_text('abC ' + str(abC))
        v.get_label_by_id('101').set_text('AbC ' + str(AbC))
        v.get_label_by_id('011').set_text('aBC ' + str(aBC))
        v.get_label_by_id('111').set_text('ABC ' + str(ABC))
        
        # Subset colors
        v.get_patch_by_id('100').set_color('c')
        v.get_patch_by_id('010').set_color('#993333')
        v.get_patch_by_id('110').set_color('blue')
        
        # Subset alphas
        v.get_patch_by_id('101').set_alpha(1.0)
        v.get_patch_by_id('011').set_alpha(1.0)
        v.get_patch_by_id('111').set_alpha(1.0)
        
        # Border styles
        c = venn3_circles(subsets=s, linestyle='solid')
        c[0].set_ls('dotted')  # Line style
        c[1].set_ls('dashed')
        c[2].set_lw(1.0)  # Line width
        plt.savefig(output_filename)
    
    if float(len(df.keys()) - 2) / 2 == 2:
        AB = 0
        for i in range(0, len(df_dict[df_dict.keys()[0]])):
            if df_dict[df_dict.keys()[0]][i] == df_dict[df_dict.keys()[1]][i]:
                AB += 1
        
        Ab = len(df_dict[df_dict.keys()[0]]) - AB
        aB = len(df_dict[df_dict.keys()[1]]) - AB
        
        s = (
            Ab,  # Ab
            aB,  # aB
            AB,  # AB
        )
        
        v = venn2(subsets=s, set_labels=('A', 'B'))
        
        # Subset labels
        v.get_label_by_id('10').set_text('Ab =' + str(Ab))
        v.get_label_by_id('01').set_text('bA =' + str(aB))
        v.get_label_by_id('11').set_text('AB =' + str(AB))
        
        # Subset colors
        v.get_patch_by_id('10').set_color('c')
        v.get_patch_by_id('01').set_color('#993333')
        v.get_patch_by_id('11').set_color('blue')
        
        # Subset alphas
        v.get_patch_by_id('10').set_alpha(0.4)
        v.get_patch_by_id('01').set_alpha(1.0)
        v.get_patch_by_id('11').set_alpha(0.7)
        
        # Border styles
        c = venn2_circles(subsets=s, linestyle='solid')
        c[0].set_ls('dashed')  # Line style
        c[0].set_lw(2.0)  # Line width
        
        plt.savefig(output_filename)


def show_help():
    print("Usage:\n")
    print("ipython venn_diagram_plots.py -c input csv_file_name -o name_of_output_png_file")


def main(argv):
    if argv:
        try:
            opts, argument = getopt.getopt(argv, "hc:o::", ["cfile=", "ofile="])
        except getopt.GetoptError:
            print('\nERROR:\n    Check your arguments', '\n')
            show_help()
            sys.exit(2)
        
        cf = ''
        of = ''
        for opt, arg in opts:
            if opt == '-h':
                show_help()
                sys.exit()
            elif opt in ("-c", "--cfile"):
                cf = arg
            elif opt in ("-o", "--ofile"):
                of = arg
        
        print('input csv filename :', cf)
        print('output png filename:', of)
        return [cf, of]


if __name__ == "__main__":
    args = main(sys.argv[1:])
    csv_file_name = args[0]
    out_name = args[1].replace('.png', '') + '.png'
    print(os.getcwd())
    read_and_plot_venn_diagram(cname=csv_file_name, output_filename=out_name)
