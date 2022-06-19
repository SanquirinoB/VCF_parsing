import ctypes

# Last stable checked version f154f8f60e1b34327eb65374e13b2ee77410f647


class Phrase(ctypes.Structure):
    _fields_ = [("m_indv", ctypes.c_ushort),
                ("m_chrom", ctypes.c_ushort),
                ("m_alele", ctypes.c_ushort),
                ("m_pos", ctypes.c_uint),
                ("m_pos_e", ctypes.c_uint),
                ("m_edit", ctypes.c_ushort),
                ("m_len", ctypes.c_uint),
                ("m_len_e", ctypes.c_uint)]


class MetaInfo(ctypes.Structure):
    _fields_ = [("m_nPhrases", ctypes.c_uint)]  # Por ahora, puede ser mas
    #    ("m_ID", ctypes.c_char * 30),
    #    ("m_relPosRef", ctypes.c_uint)]


class MetaRef(ctypes.Structure):
    _fields_ = [("m_ID", ctypes.c_uint),
               ("m_nBases", ctypes.c_uint),
               ("m_relPos", ctypes.c_uint)]