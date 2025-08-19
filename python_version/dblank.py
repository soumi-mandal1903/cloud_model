class StringUtils:
    @staticmethod
    def dblank(s):
        """
        Returns the number of characters up to and including the last non-blank character.
        Equivalent to Fortran's dblank.
        """
        ns = len(s)
        while ns > 1 and s[ns-1] == ' ':
            ns -= 1
        return ns

# Example usage:
# s = "hello   "
# ns = StringUtils.dblank(s)
# print(ns)  # Output: 5