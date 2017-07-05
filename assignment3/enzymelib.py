
class Enzyme:
    ###########################################################################
    # ENZYMES
    ###########################################################################
    # Trypsin
    def trypsin(self, seq):
        splitmatrix = []
        seq_piece = ""
        for i in range(len(seq)):
            seq_piece += seq[i]
            if (seq[i] == 'K' or seq[i] == 'R'):
                if i < len(seq) - 1 and seq[i + 1] == 'P':
                    continue
                else:
                    splitmatrix.append(seq_piece)
                    seq_piece = ""
        return splitmatrix

    ###########################################################################
    # ENZYME REGISTRATION
    ###########################################################################
    __functions = {'trypsin': trypsin}

    ###########################################################################
    # ENZYME INIT
    ###########################################################################
    def __init__(self, enz_ident):
        if enz_ident in self.__functions:
            self.name = enz_ident
            self.fnc = self.__functions[enz_ident]
        else:
            print('Invalid digestion enzyme! '
                  'Using default enzyme "trypsin" (Trypsin)')
            self.name = 'trypsin'
            self.fnc = self.__functions['trypsin']

    ###########################################################################
    # ENZYME FUNCTIONS
    ###########################################################################
    def cleave(self, seq):
        return self.fnc(self, seq)
