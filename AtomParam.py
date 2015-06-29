""" Stores the relevant parameters for the used atoms """

""" Self energies of the atoms """
e0 = {'C':   0.0,
      'B':   1.7,           #verified
      'N':  -1.9,           #verified
      'O':   2.6,           #verified
      'S':  -2.2}  
       
""" Neirest neighbour interaction between atoms """
tm = {'CC':-2.8, 'CN': -2.8, 'NC': -2.8,    #verified
                 'CB': -2.8, 'BC': -2.8,    #verified
                 'CO': -2.8, 'OC': -2.8,    #verifed
                 'CS':  6.6, 'SC':  6.6}    # To be verified

""" Colors for displaying the atoms """
cAtom = {'C': 'k',
         'B': 'gray',
         'N': 'g',
         'O': 'b',
         'S': 'y'}

""" Initial interaction values """
gLeft  = 0.2
gRight = 0.2

""" Hubbard parameter """
Hubbard = 8
