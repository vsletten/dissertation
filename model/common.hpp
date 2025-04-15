/*****************************************************************************
 * Here is the site classification scheme:
 *
 *   EDGE   used to identify site at edge of lattice
 *
 *   Al sites
 *   --------
 *   100   empty
 *   101   Al(OH,H2O)0
 *   102   Al(OH,H2O)1
 *   103   Al(OH,H2O)2
 *   104   Al(OH,H2O)3
 *   105   Al(OH,H2O)4
 *   106   Al(OH,H2O)5
 *   107   Al(OH,H20)6
 *   199   Si(OH)4
 *
 *   Si sites
 *   --------
 *   200   empty
 *   201   Si(OH)0
 *   202   Si(OH)1
 *   203   Si(OH)2
 *   204   Si(OH)3
 *   205   Si(OH)4
 *   299   Al(OH,H2O)6
 *
 *   Si-O-Si type O sites
 *   --------------------
 *   300    empty
 *   301    Si-O-Si
 *   302    Si-OH HO-Si
 *   303    Si-OH
 *
 *   Si-O<Al2 type O sites
 *   ---------------------
 *   400    empty
 *   401    Si-O<Al2
 *   402    Si-OH HO<Al2
 *   403    Si-OH HO-Al H2O-Al
 *   404          HO<Al2
 *   405          HO-Al H2O-Al
 *   406    Si-OH-Al
 *   407    Si-OH HO-Al
 *   408    Si-OH
 *   409          HO-Al
 *   410    Si-OH-Al HO-Al
 *
 *   Al-OH-Al type O(H) sites
 *   ------------------------
 *   500    empty
 *   501    Al-OH-Al
 *   502    Al-OH H2O-Al
 *   503    Al-OH
 *
 *****************************************************************************
 * Within this scheme, there are 16 possible reactions (not
 * including adsorption/desorption), with some duplicated effort:
 *
 *  R0. 301 --> 302   Si-O-Si + H2O = Si-OH + Si-OH
 *  R1. 302 --> 301
 *  R2. 401 --> 402   Si-O-Al2 + H2O = Si-OH + Al-OH-Al
 *  R3. 402 --> 401
 *  R4. 401 --> 410   Si-O-Al2 + H2O = Si-OH-Al + Al-OH  
 *  R5. 410 --> 401
 *  R6. 402 --> 403   Si-OH + Al-OH-Al + H2O = Si-OH + Al-OH + Al-H2O  
 *  R7. 403 --> 402
 *  R8. 410 --> 403   Si-OH-Al + Al-OH + H2O = Si-OH + Al-OH + Al-H2O  
 *  R9. 403 --> 410
 * R10. 406 --> 407   Si-OH-Al + H2O = Si-OH + Al-H2O  
 * R11. 407 --> 406
 * R12. 404 --> 405   Al-OH-Al + H2O = Al-OH + Al-H2O
 * R13. 405 --> 404  
 * R14. 501 --> 502   Al-OH-Al + H2O = Al-OH + Al-H2O
 * R15. 502 --> 501
 *
 * 
 * Transitions possible upon adsorption:
 * 
 *   1. 100 --> 107  (Al)
 *   2. 200 --> 205  (Si)
 *   3. 300 --> 303  (Si)
 *   4. 303 --> 302  (Si)
 *   5. 400 --> 408  (Si)
 *   6. 400 --> 409  (Al)
 *   7. 404 --> 402  (Si)
 *   8. 405 --> 403  (Si)
 *   9. 406 --> 410  (Al)
 *  10. 407 --> 403  (Al)
 *  11. 408 --> 407  (Al)
 *  12. 409 --> 405  (Al)
 *  13. 409 --> 407  (Si)
 *  14. 500 --> 503  (Al)
 *  15. 503 --> 502  (Al)
 *
 * Transitions possible upon desorption:
 * 
 *   1. 107 --> 100  (Al)
 *   2. 205 --> 200  (Si)
 *   3. 303 --> 300  (Si)
 *   4. 302 --> 303  (Si)
 *   5. 402 --> 404  (Si)
 *   6. 403 --> 405  (Si)
 *   7. 410 --> 406  (Al)
 *   8. 403 --> 407  (Al)
 *   9. 407 --> 408  (Al)
 *  10. 405 --> 409  (Al)
 *  11. 407 --> 409  (Si)
 *  12. 409 --> 400  (Al)
 *  13. 408 --> 400  (Si)
 *  14. 502 --> 503  (Al)
 *  15. 503 --> 500  (Al)
 *
 * Next-neighbor state classification scheme (index, x, y):
 *
 *  Adsorption:
 *
 *   100's:  x = (double 500 occupied) ? 1 : 0
 *           y = number of occupied 400's
 *               (0 - 2)
 *       index = x * 3 + y
 *
 *   200's:  x = (400 occupied) ? 1 : 0
 *           y = number of occupied 300's
 *               (0 - 3)
 *       index = x * 4 + y
 *
 *  Desorption:
 *
 *   100's:  x = (double 500 occupied) ? 1 : 0
 *           y = number of occupied 400's
 *               (0 - 2)
 *       index = x * 3 + y
 *
 *   200's:  x = (400 occupied) ? 1 : 0
 *           y = number of occupied 300's
 *               (0 - 3)
 *       index = x * 4 + y
 * 
 *  Hydrolysis (and reverse):
 *
 *   300's:  x = number of (402, 403, 406, 407, 408, 410) on Si
 *               (0 - 2)
 *           y = number of (302, 303) on Si
 *               (0 - 4)
 *       index = x * 5 + y
 *
 *   400's:  x = number of (302, 303) on Si
 *               (0 - 3)
 *           y = number of (403, 405, 406, 407, 409, 410, 502, 503) on Al
 *               (0 - 9)
 *       index = x * 10 + y
 *
 *   500's:  x = (other bridge hydrolyzed) ? 1 : 0
 *           y = number of (403, 405, 406, 407, 409, 410, 502, 503) on Al
 *               (0 - 8)
 *       index = x * 9 + y
 *
 *
 *****************************************************************************/

#ifndef common_h
#define common_h

#ifndef FALSE
  #define FALSE 0
  #define TRUE !FALSE
#endif

typedef enum {UNREACHABLE, ENQUEUED, DISCOVERED} Color;



#endif

