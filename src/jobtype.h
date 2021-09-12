/* Internal handles for job types */
# define J_HF   1
# define J_GW   2
# define J_RPA  3
# define J_HFE  4
# define J_GT   5
# define J_SUS  6
# define J_SUSR 7
# define J_DIEL 8
# define J_SCR  9
# define J_SCRW 10
# define J_GOLD 11
# define J_SX   12
# define J_COSX 13
# define J_PBE0 14
# define J_KS   15

/* Job types for spectra */
# define J_SPEC [ J_SUS,J_SUSR,J_DIEL,J_SCR ]

/* Job types that require bare exchange */
# define J_EXCH [ J_HF,J_GW,J_RPA,J_HFE,J_PBE0 ]
