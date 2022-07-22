import ctypes


class Phrase(ctypes.Structure):
    _fields_ = [("m_indv", ctypes.c_longlong),
                ("m_chrom", ctypes.c_longlong),
                ("m_alele", ctypes.c_longlong),
                ("m_pos", ctypes.c_longlong),
                ("m_pos_e", ctypes.c_longlong),
                ("m_edit", ctypes.c_longlong),
                ("m_len", ctypes.c_longlong),
                ("m_len_e", ctypes.c_longlong)]


class MetaInfo(ctypes.Structure):
    _fields_ = [("m_nPhrases", ctypes.c_longlong)]  # Por ahora, puede ser mas
    #    ("m_ID", ctypes.c_char * 30),
    #    ("m_relPosRef", ctypes.c_longlong)]


class MetaRef(ctypes.Structure):
    _fields_ = [("m_ID", ctypes.c_longlong),
               ("m_nBases", ctypes.c_longlong),
               ("m_relPos", ctypes.c_longlong)]