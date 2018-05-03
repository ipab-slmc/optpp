
/*********************************************************************
 *
 * Common block structures for constraints and parallel configuration.
 *
 *********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

struct conbcmni {
    int ncb;
    int ncni;
};

struct pdscon {
    int me;
    int nproc;
};

#ifdef __cplusplus
}
#endif

