import ctypes

# Last stable checked version f154f8f60e1b34327eb65374e13b2ee77410f647


class Phrase(ctypes.Structure):
    _fields_ = [("m_indv", ctypes.c_ushort),
               ("m_chrom", ctypes.c_ubyte),
               ("m_alele", ctypes.c_ubyte),
               ("m_pos", ctypes.c_uint),
               ("m_pos_e", ctypes.c_uint),
               ("m_edit", ctypes.c_ushort),
               ("m_len", ctypes.c_uint),
               ("m_len_e", ctypes.c_uint)]

if __name__ == "__main__":
    print(ctypes.sizeof(Phrase))