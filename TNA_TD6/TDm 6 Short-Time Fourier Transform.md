> **Objective**: Verify COLA property, OLA and WOLA method

## A-/ Check constant-overlap-add criteria of following examples
- Rectangular window at 0% overlap (hop size R = window size M )
- Bartlett window at 50% overlap ( R ~= M/2 ) (Since normally M is odd, "R ~=
M/2" means "R=(M-1)/2", etc.)
- Hamming window at 50% overlap ( R ~= M/2 )
- Rectangular window at 50% overlap ( R ~= M/2 )
- Hamming window at 75% overlap ( R = M/4 = 25% hop size)
- Blackman family at 2/3 overlap (1/3 hop size); e.g.,
blackman(33,'periodic')
## B-/ Implement STFT method (OLA) using hamming wodow w/50% Ovelap and verify the reconstruction in temporal domain
- x(n) = somme(xk(n))
- OLA principle including STFT and ISTFT.
## C-/ Implement WOLA method and propose a unique widow for both analysis and synthesis