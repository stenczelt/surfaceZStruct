class CoordinationNumbers:
    """
    class to store maximum allowed coordination number for each atom.
    max coordination number is retrieved based on atomic numbers
    """

    # def __init__(self):

    def getMaxCoordNum(self, atomicNumber):
        atomicNum_coordNum_dictionary = {
            1: 1,  # H
            2: 0,  # He
            3: 1,  # Li
            5: 3,  # B
            6: 4,  # C
            7: 3,  # N
            8: 2,  # O
            9: 1,  # F
            22: 6,  # Ti
        }
        return atomicNum_coordNum_dictionary[atomicNumber]
