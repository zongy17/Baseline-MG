In the following context, one rankfile represents one test of a specific MPI/OpenMP configuration. 
Multigrids of Struct and SStruct interface need 3D process partitions, and thus each row in their tables is a test.
Since unstructured multigrid (i.e., BoomerAMG) only needs 1D process partitions, multiple tests (i.e., multiple rankfiles) are written together in a same row to save context space.
# Problem: rhd
case_name = "LASER"
| hypre-MG | SMG & PFMG | BoomerAMG |
|---|:--------:|:---:|
| executable | base_struct_f64.exe | base_IJ_f64.exe |
## Case: size of 128x128x128
This case is for strong scalability tests.

cx = 128, cy = 128, cz = 128
### SMG & PFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1 | N0.5P64| 4 | 4 | 4 |
| 1 | N1P64  | 4 | 4 | 4 |
| 1 | N1P128 | 4 | 8 | 4 |
| 2 | N2P64  | 4 | 4 | 4 |
| 2 | N2P128 | 4 | 8 | 4 |
| 2 | N2P256 | 8 | 8 | 4 |
| 4 | N4P64  | 4 | 4 | 4 |
| 4 | N4P128 | 4 | 8 | 4 |
| 4 | N4P256 | 8 | 8 | 4 |
| 4 | N4P512 | 8 | 8 | 8 |
| 8 | N8P64  | 4 | 4 | 4 |
| 8 | N8P128 | 4 | 8 | 4 |
| 8 | N8P256 | 8 | 8 | 4 |
| 8 | N8P512 | 8 | 8 | 8 |
| 8 | N8P1024| 8 |16 | 8 |
|16 |N16P256 | 8 | 8 | 4 |
|16 |N16P512 | 8 | 8 | 8 |
|16 |N16P1024| 8 |16 | 8 |
|16 |N16P2048|16 |16 | 8 |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1 | N0.5P4, N0.5P8, N0.5P16, N0.5P32, N0.5P64|
| 1 | N1P8  , N1P16 , N1P32  , N1P64  , N1P128 |
| 2 | N2P16 , N2P32 , N2P64  , N2P128 , N2P256 |
| 4 | N4P32 , N4P64 , N4P128 , N4P256 , N4P512 |
| 8 | N8P64 , N8P128, N8P256 , N8P512 , N8P1024|
|16 | N16P128, N16P256, N16P512 ,N16P1024, N16P2048|

# Problem: oil
case_name = "OIL-IMPEC"
| hypre-MG | BoomerAMG |
|---|:--------:|
| executable | base_IJ_f64.exe |
## Case: size of 640x64x768
This case is for strong scalability tests.

cx = 640, cy = 64, cz = 768
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 2 | N2P24 , N2P40 , N2P48  , N2P80 , N2P120, N2P240 |
| 4 | N4P48 , N4P80 , N4P96  , N4P160, N4P240, N4P480 |
| 8 | N8P96 , N8P160, N8P192 , N8P320, N8P480, N8P960 |
|16 | N16P192, N16P320, N16P384 , N16P640, N16P960, N16P1920 |
|32 | N32P384, N32P640, N32P768 , N32P1280, N32P1920, N32P3840 |


# Problem: weather
case_name = "GRAPES"
| hypre-MG | SMG & PFMG | BoomerAMG |
|---|:--------:|:---:|
| executable | base_struct_f32.exe | base_IJ_f32.exe / base_IJ_f32_i64.exe |

**Notice that base_IJ_f32.exe is for Case 1, 2, 3, while base_IJ_f32_i64.exe for Case 4 whose nonzero indices of CSR exceed $2^{31}-1$.**
## Case 1: size of 480x288x72
This case is for weak scalability tests.

cx = 480, cy = 288, cz = 72
### SMG & PFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1    |  N1P120  | 10 | 12 | 1  |
| 1    |  N1P60   | 10 |  6 | 1  |
| 1    |  N1P40   |  5 |  8 | 1  |
| 1    |  N1P24   |  4 |  6 | 1  |
| 1    |  N1P20   |  5 |  4 | 1  |
| 2    |  N2P240  | 20 | 12 | 1  |
| 2    |  N2P120  | 10 | 12 | 1  |
| 2    |  N2P80   | 10 |  8 | 1  |
| 2    |  N2P48   |  8 |  6 | 1  |
| 2    |  N2P40   |  5 |  8 | 1  |
| 4    |  N4P480  | 20 | 24	| 1  |
| 4    |  N4P240  | 20 | 12 | 1  |
| 4    |  N4P160  | 10 | 16 | 1  |
| 4    |  N4P96   |  8 | 12 | 1  |
| 4    |  N4P80   | 10 |  8 | 1  |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1    |  N1P120, N1P60 , N1P40, N1P24, N1P20 |
| 2    |  N2P240, N2P120, N2P80, N2P48, N2P40 |
| 4    |  N4P480, N4P240, N4P160, N4P96, N4P80 |

## Case 2: size of 960x576x72
This case is for weak scalability tests.

cx = 960, cy = 576, cz = 72
### SMG & PFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1    |  N1P120  | 10 | 12 | 1  |
| 1    |  N1P60   | 10 |  6 | 1  |
| 1    |  N1P40   |  5 |  8 | 1  |
| 1    |  N1P24   |  4 |  6 | 1  |
| 1    |  N1P20   |  5 |  4 | 1  |
| 2    |  N2P240  | 20 | 12 | 1  |
| 2    |  N2P120  | 10 | 12 | 1  |
| 2    |  N2P80   | 10 |  8 | 1  |
| 2    |  N2P48   |  8 |  6 | 1  |
| 2    |  N2P40   |  5 |  8 | 1  |
| 4    |  N4P480  | 20 | 24	| 1  |
| 4    |  N4P240  | 20 | 12 | 1  |
| 4    |  N4P160  | 10 | 16 | 1  |
| 4    |  N4P96   |  8 | 12 | 1  |
| 4    |  N4P80   | 10 |  8 | 1  |
| 8    |  N8P960  | 40 | 24 | 1  |
| 8    |  N8P480  | 20 | 24 | 1  |
| 8    |  N8P320  | 20 | 16 | 1  |
| 8    |  N8P192  | 16 | 12 | 1  |
| 8    |  N8P160  | 10 | 16 | 1  |
| 16   | N16P1920 | 40 | 48 | 1  |
| 16   | N16P960  | 40 | 24 | 1  |
| 16   | N16P640  | 20 | 32 | 1  |
| 16   | N16P384  | 16 | 24 | 1  |
| 16   | N16P320  | 20 | 16 | 1  |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1    |  N1P120, N1P60 , N1P40, N1P24, N1P20 |
| 2    |  N2P240, N2P120, N2P80, N2P48, N2P40 |
| 4    |  N4P480, N4P240, N4P160, N4P96, N4P80 |
| 8    |  N8P960, N8P480, N8P320, N8P192, N8P160 |
| 16   | N16P1920, N16P960, N16P640, N16P384, N16P320 |

## Case 3: size of 1920x1152x72
This case is for weak scalability tests.

cx = 1920, cy = 1152, cz = 72
### SMG & PFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1    |  N1P120  | 10 | 12 | 1  |
| 1    |  N1P60   | 10 |  6 | 1  |
| 1    |  N1P40   |  5 |  8 | 1  |
| 1    |  N1P24   |  4 |  6 | 1  |
| 1    |  N1P20   |  5 |  4 | 1  |
| 2    |  N2P240  | 20 | 12 | 1  |
| 2    |  N2P120  | 10 | 12 | 1  |
| 2    |  N2P80   | 10 |  8 | 1  |
| 2    |  N2P48   |  8 |  6 | 1  |
| 2    |  N2P40   |  5 |  8 | 1  |
| 4    |  N4P480  | 20 | 24	| 1  |
| 4    |  N4P240  | 20 | 12 | 1  |
| 4    |  N4P160  | 10 | 16 | 1  |
| 4    |  N4P96   |  8 | 12 | 1  |
| 4    |  N4P80   | 10 |  8 | 1  |
| 8    |  N8P960  | 40 | 24 | 1  |
| 8    |  N8P480  | 20 | 24 | 1  |
| 8    |  N8P320  | 20 | 16 | 1  |
| 8    |  N8P192  | 16 | 12 | 1  |
| 8    |  N8P160  | 10 | 16 | 1  |
| 16   | N16P1920 | 40 | 48 | 1  |
| 16   | N16P960  | 40 | 24 | 1  |
| 16   | N16P640  | 20 | 32 | 1  |
| 16   | N16P384  | 16 | 24 | 1  |
| 16   | N16P320  | 20 | 16 | 1  |
| 32   | N32P3840 | 80 | 48 | 1  |
| 32   | N32P1920 | 40 | 48 | 1  |
| 32   | N32P1280 | 40 | 32 | 1  |
| 32   | N32P768  | 32 | 24 | 1  |
| 32   | N32P640  | 20 | 32 | 1  |
| 32   | N32P384  | 16 | 24 | 1  |
| 64   | N64P7680 | 80 | 96 | 1  |
| 64   | N64P3840 | 80 | 48 | 1  |
| 64   | N64P2560 | 40 | 64 | 1  |
| 64   | N64P1536 | 32 | 48 | 1  |
| 64   | N64P1280 | 40 | 32 | 1  |
| 64   | N64P768  | 32 | 24 | 1  |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1    |  N1P120, N1P60 , N1P40, N1P24, N1P20 |
| 2    |  N2P240, N2P120, N2P80, N2P48, N2P40 |
| 4    |  N4P480, N4P240, N4P160, N4P96, N4P80 |
| 8    |  N8P960, N8P480, N8P320, N8P192, N8P160 |
| 16   | N16P1920, N16P960, N16P640, N16P384, N16P320 |
| 32   | N32P3840, N32P1920, N32P1280, N32P768, N32P640, N32P384 |
| 64   | N64P7680, N64P3840, N64P2560, N64P1536, N64P1280, N64P768 |

## Case 4: size of 3840x2304x72
This case is for weak and strong scalability tests.

cx = 3840, cy = 2304, cz = 72
### SMG & PFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 2    |  N2P240  | 20 | 12 | 1  |
| 2    |  N2P120  | 10 | 12 | 1  |
| 2    |  N2P80   | 10 |  8 | 1  |
| 2    |  N2P48   |  8 |  6 | 1  |
| 2    |  N2P40   |  5 |  8 | 1  |
| 4    |  N4P480  | 20 | 24	| 1  |
| 4    |  N4P240  | 20 | 12 | 1  |
| 4    |  N4P160  | 10 | 16 | 1  |
| 4    |  N4P96   |  8 | 12 | 1  |
| 4    |  N4P80   | 10 |  8 | 1  |
| 8    |  N8P960  | 40 | 24 | 1  |
| 8    |  N8P480  | 20 | 24 | 1  |
| 8    |  N8P320  | 20 | 16 | 1  |
| 8    |  N8P192  | 16 | 12 | 1  |
| 8    |  N8P160  | 10 | 16 | 1  |
| 16   | N16P1920 | 40 | 48 | 1  |
| 16   | N16P960  | 40 | 24 | 1  |
| 16   | N16P640  | 20 | 32 | 1  |
| 16   | N16P384  | 16 | 24 | 1  |
| 16   | N16P320  | 20 | 16 | 1  |
| 32   | N32P3840 | 80 | 48 | 1  |
| 32   | N32P1920 | 40 | 48 | 1  |
| 32   | N32P1280 | 40 | 32 | 1  |
| 32   | N32P768  | 32 | 24 | 1  |
| 32   | N32P640  | 20 | 32 | 1  |
| 32   | N32P384  | 16 | 24 | 1  |
| 64   | N64P7680 | 80 | 96 | 1  |
| 64   | N64P3840 | 80 | 48 | 1  |
| 64   | N64P2560 | 40 | 64 | 1  |
| 64   | N64P1536 | 32 | 48 | 1  |
| 64   | N64P1280 | 40 | 32 | 1  |
| 64   | N64P768  | 32 | 24 | 1  |
| 128  |N128P15360| 160| 96 | 1  |
| 128  |N128P7680 | 80 | 96 | 1  |
| 128  |N128P5120 | 80 | 64 | 1  |
| 128  |N128P3072 | 64 | 48 | 1  |
| 128  |N128P2560 | 40 | 64 | 1  |
| 128  |N128P1536 | 32 | 48 | 1  |
| 256  |N256P30720| 160|192 | 1  |
| 256  |N256P15360| 160| 96 | 1  |
| 256  |N256P10240| 80 |128 | 1  |
| 256  |N256P6144 | 64 | 96 | 1  |
| 256  |N256P5120 | 80 | 64 | 1  |
| 256  |N256P3072 | 64 | 48 | 1  |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 2    |  N2P240, N2P120, N2P80, N2P48, N2P40 |
| 4    |  N4P480, N4P240, N4P160, N4P96, N4P80 |
| 8    |  N8P960, N8P480, N8P320, N8P192, N8P160 |
| 16   | N16P1920, N16P960, N16P640, N16P384, N16P320 |
| 32   | N32P3840, N32P1920, N32P1280, N32P768, N32P640, N32P384 |
| 64   | N64P7680, N64P3840, N64P2560, N64P1536, N64P1280, N64P768 |
| 128  |N128P15360, N128P7680, N128P5120, N128P3072, N128P2560, N128P1536 |
| 256  |N256P30720, N256P15360, N256P10240, N256P6144, N256P5120, N256P3072 |

# Problem: rhd-3T
case_name = "LASER-3T"
| hypre-MG | SysPFMG | BoomerAMG |
|---|:--------:|:---:|
| executable | base_sstruct_f64.exe | base_IJ_f64.exe |
## Case: size of 128x128x128
This case is for strong scalability tests.

cx = 128, cy = 128, cz = 128
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1 | N0.5P64| 4 | 4 | 4 |
| 1 | N1P64  | 4 | 4 | 4 |
| 1 | N1P128 | 4 | 8 | 4 |
| 2 | N2P64  | 4 | 4 | 4 |
| 2 | N2P128 | 4 | 8 | 4 |
| 2 | N2P256 | 8 | 8 | 4 |
| 4 | N4P64  | 4 | 4 | 4 |
| 4 | N4P128 | 4 | 8 | 4 |
| 4 | N4P256 | 8 | 8 | 4 |
| 4 | N4P512 | 8 | 8 | 8 |
| 8 | N8P64  | 4 | 4 | 4 |
| 8 | N8P128 | 4 | 8 | 4 |
| 8 | N8P256 | 8 | 8 | 4 |
| 8 | N8P512 | 8 | 8 | 8 |
| 8 | N8P1024| 8 |16 | 8 |
|16 |N16P256 | 8 | 8 | 4 |
|16 |N16P512 | 8 | 8 | 8 |
|16 |N16P1024| 8 |16 | 8 |
|16 |N16P2048|16 |16 | 8 |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1 | N0.5P4, N0.5P8, N0.5P16, N0.5P32, N0.5P64|
| 1 | N1P8  , N1P16 , N1P32  , N1P64  , N1P128 |
| 2 | N2P16 , N2P32 , N2P64  , N2P128 , N2P256 |
| 4 | N4P32 , N4P64 , N4P128 , N4P256 , N4P512 |
| 8 | N8P64 , N8P128, N8P256 , N8P512 , N8P1024|
|16 | N16P128, N16P256, N16P512 ,N16P1024, N16P2048|

# Problem: oil-4C
case_name = "OIL-FIM"
| hypre-MG | BoomerAMG |
|---|:--------:|
| executable | base_IJ_f64.exe |
## Case: size of 320x64x384
This case is for strong scalability tests.

cx = 320, cy = 64, cz = 384
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1    |  N1P120, N1P60 , N1P40, N1P24, N1P20 |
| 2    |  N2P240, N2P120, N2P80, N2P48, N2P40 |
| 4    |  N4P480, N4P240, N4P160, N4P96, N4P80 |
| 8    |  N8P960, N8P480, N8P320, N8P192, N8P160 |
| 16   | N16P1920, N16P960, N16P640, N16P384, N16P320 |
| 32   | N32P3840, N32P1920, N32P1280, N32P768, N32P640, N32P384 |
| 64   | N64P7680, N64P3840, N64P2560, N64P1536, N64P1280, N64P768 |

# Problem: solid-3D
case_name = "SOLID"
| hypre-MG | SysPFMG | BoomerAMG |
|---|:--------:|:---:|
| executable | base_sstruct_f64.exe | base_IJ_f64.exe |
## Case 1: 80x96x64
This case is for weak scalability tests.

cx = 80, cy = 96, cz = 64
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1 | N0.5P30| 5 | 3 | 2 |
| 1 | N0.5P60| 5 | 3 | 4 |
| 1 | N1P30  | 5 | 3 | 2 |
| 1 | N1P60  | 5 | 3 | 4 |
| 1 | N1P120 | 5 | 3 | 8 |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1 | N0.5P60, N0.5P30, N0.5P20, N0.5P10, N0.5P6 |
| 1 | N1P120, N1P60 , N1P40, N1P24, N1P20, N1P12 |

## Case 2: 80x96x128
This case is for weak scalability tests.

cx = 80, cy = 96, cz = 128
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1 | N1P60  | 5 | 6 | 2 |
| 1 | N1P120 | 5 | 6 | 4 |
| 2 | N2P240 | 5 | 6 | 8 |
| 2 | N2P120 | 5 | 6 | 4 |
| 2 | N2P60  | 5 | 3 | 4 |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 1 | N1P120, N1P60 , N1P40, N1P24, N1P20, N1P12 |
| 2 | N2P240, N2P120, N2P80, N2P48, N2P40, N2P24 |

## Case 3: 160x96x128
This case is for weak scalability tests.

cx = 160, cy = 96, cz = 128
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 2 | N2P240 | 5 | 6 | 8 |
| 2 | N2P120 | 5 | 6 | 4 |
| 2 | N2P60  | 5 | 3 | 4 |
| 2 | N2P60  | 5 | 6 | 2 |
| 4 | N4P480 |10 | 6 | 8 |
| 4 | N4P240 | 5 | 6 | 8 |
| 4 | N4P120 | 5 | 6 | 4 |
### BoomerAMG
For unstructured AMG, only specifying rankfiles is enough, because only 1D process partition is needed.
| Node | Rankfiles |
|---|:--------:|
| 2 | N2P240, N2P120, N2P80, N2P48, N2P40, N2P24 |
| 4 | N4P480, N4P240, N4P160, N4P96, N4P80, N4P48 |

## Case 4: 160x192x128
This case is for weak and strong scalability tests.

cx = 160, cy = 192, cz = 128
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 1 | N1P120 | 5 | 6 | 4 |
| 1 | N1P60 | 5 | 3 | 4 |
| 1 | N1P30 | 5 | 3 | 2 |
| 2 | N2P240| 5 | 6 | 8 |
| 2 | N2P120| 5 | 6 | 4 |
| 2 | N2P60 | 5 | 3 | 4 |
| 4 | N4P480|10 | 6 | 8 |
| 4 | N4P240| 5 | 6 | 8 |
| 4 | N4P120| 5 | 6 | 4 |
| 8 | N8P960|10 |12 | 8 |
| 8 | N8P480|10 | 6 | 8 |
| 8 | N8P240| 5 | 6 | 8 |
|16 |N16P960|10 |12 | 8 |
|16 |N16P480|10 | 6 | 8 |
|16 |N16P480|10 |12 | 4 |
|32 |N32P3840|10|24 |16 |
|32 |N32P1920|10|12 |16 |
|32 |N32P1280|10| 8 |16 |
|32 |N32P960 |10|12 | 8 |
|32 |N32P640 |10| 8 | 8 |
### BoomerAMG
| Node | Rankfiles |
|---|:--------:|
| 1 | N1P120, N1P60 , N1P40, N1P24, N1P20, N1P12 |
| 2 | N2P240, N2P120, N2P80, N2P48, N2P40, N2P24 |
| 4 | N4P480, N4P240, N4P160, N4P96, N4P80, N4P48 |
| 8 | N8P960, N8P480, N8P320, N8P192, N8P160, N8P96 |
|16 | N16P1920, N16P960, N16P640, N16P384, N16P320, N16P192 |
|30 | N30P3600, N30P1800, N30P1200, N30P720, N30P600, N30P360 |
|32 | N32P3840, N32P1920, N32P1280, N32P768, N32P640, N32P384 |

## Case 5: 160x192x256
This case is for weak scalability tests.

cx = 160, cy = 192, cz = 256
### SysPFMG
| Node | Rankfile | px | py | pz |
|---|:--------:|---:|----|----|
| 8 | N8P960  | 10 | 12 | 8 |
| 8 | N8P480  | 10 | 6  | 8 |
| 8 | N8P240  | 5  | 6  | 8 |
| 8 | N8P120  | 5  | 6  | 4 |
|16 |N16P1920 | 10 | 12 |16 |
|16 |N16P960  | 10 | 12 | 8 |
|16 |N16P480  | 10 | 6  | 8 |
|16 |N16P240  | 5  | 6  | 8 |
### BoomerAMG
| Node | Rankfiles |
|---|:--------:|
| 8 | N8P960, N8P480, N8P320, N8P192, N8P160, N8P96 |
|16 | N16P1920, N16P960, N16P640, N16P384, N16P320, N16P192 |