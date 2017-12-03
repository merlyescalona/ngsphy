class NGSphyException(Exception):
    """
    Exception raised for errors of the NGSphy program.
    ----------------------------------------------------------------------------
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message, time):
        self.expression = expression
        self.message = message
        self.time= time

class NGSphyExitException(Exception):
    """
    Exception raised for ending of the NGSphy process.
    ----------------------------------------------------------------------------
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, expression, message, time):
        self.expression = expression
        self.message = message
        self.time= time
