#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//FFTReal
#include <cassert>
#include <cmath>

//rose_image
#include "TImage.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TROOT.h"
#include "TColor.h"
#include "TAttImage.h"
#include "TEnv.h"

TGraph* gMutagenTestGraph = NULL;
double* gMutagenOrdinates = NULL;
//bool gMutagenGraphOn = true;

void mutagenInitGraphics(int numMaxOrdinates)
{
    if (NULL != gMutagenOrdinates) {
        delete[] gMutagenOrdinates;
        gMutagenOrdinates = NULL;
    }
    
    gMutagenOrdinates = new double[numMaxOrdinates];
    for (int i = 0; i < numMaxOrdinates; i++) {
        gMutagenOrdinates[i] = (double)i;
    }
    
    //if (NULL == pMutagenTestGraph) {
    //    delete pMutagenTestGraph;
    //    pMutagenTestGraph = NULL;
    //}
}

void mutagenSafeDelete(void* p, bool array)
{
    if (NULL != p) {
        if (array)
            delete[] p;
        else
            delete p
    }
}

struct WavInfo
{
  FILE* pHandle;
  unsigned short wNumChannels;
  unsigned int lSampleRate;
  unsigned short wBytesPerSample;
  unsigned int lNumFrames;
  unsigned int lNumDataBytes;
};

class WaveFile
{
  static int openRead(char* pszPath, WavInfo* pWavInfo);
  static FILE* openWrite(char* pszPath);
  static int read(WavInfo* pWavInfo, float* pwData);
  static int write(WavInfo* pWavInfo, float* pwData);
};

WaveFile gWaveFile;
/*
The canonical WAVE format starts with the RIFF header:

0         4   ChunkID          Contains the letters "RIFF" in ASCII form
                               (0x52494646 big-endian form).
4         4   ChunkSize        36 + SubChunk2Size, or more precisely:
                               4 + (8 + SubChunk1Size) + (8 + SubChunk2Size)
                               This is the size of the rest of the chunk 
                               following this number.  This is the size of the 
                               entire file in bytes minus 8 bytes for the
                               two fields not included in this count:
                               ChunkID and ChunkSize.
8         4   Format           Contains the letters "WAVE"
                               (0x57415645 big-endian form).

The "WAVE" format consists of two subchunks: "fmt " and "data":
The "fmt " subchunk describes the sound data's format:

12        4   Subchunk1ID      Contains the letters "fmt "
                               (0x666d7420 big-endian form).
16        4   Subchunk1Size    16 for PCM.  This is the size of the
                               rest of the Subchunk which follows this number.
20        2   AudioFormat      PCM = 1 (i.e. Linear quantization)
                               Values other than 1 indicate some 
                               form of compression.
22        2   NumChannels      Mono = 1, Stereo = 2, etc.
24        4   SampleRate       8000, 44100, etc.
28        4   ByteRate         == SampleRate * NumChannels * BitsPerSample/8
32        2   BlockAlign       == NumChannels * BitsPerSample/8
                               The number of bytes for one sample including
                               all channels. I wonder what happens when
                               this number isn't an integer?
34        2   BitsPerSample    8 bits = 8, 16 bits = 16, etc.
          2   ExtraParamSize   if PCM, then doesn't exist
          X   ExtraParams      space for extra parameters

The "data" subchunk contains the size of the data and the actual sound:

36        4   Subchunk2ID      Contains the letters "data"
                               (0x64617461 big-endian form).
40        4   Subchunk2Size    == NumSamples * NumChannels * BitsPerSample/8
                               This is the number of bytes in the data.
                               You can also think of this as the size
                               of the read of the subchunk following this 
                               number.
44        *   Data             The actual sound data.
*/

// Checks
// W Bytes 0-3 say "RIFF"
// W Subchunk1Size + Subchunk2Size + 20 == ChunkSize
// W ChunkSize + 8 == FileSize
// W Bytes 8-11 say "WAVE"
// W Bytes 12-15 say "fmt "
// W Subchunk1Size == 16
// W AudioFormat == 1
// E NumSamples % NumChannels == 0
// E BlockAlign * SampleRate == ByteRate
// E BitsPerSample % 8 == 0
// E NumChannels * BitsPerSample / 8 == BlockAlign
// W Bytes 36-39 say "data"
// E Subchunk2Size % BlockAlign == 0

int WaveFile::openRead(char* pszPath, WavInfo* pWavInfo)
{
    if (NULL == pWavInfo)
        return -1;

    pWavInfo->pHandle = fopen(pszPath, "rb");
    if (NULL == pWavInfo->pHandle)
        return -2;

    fseek(pWavInfo->pHandle, 0, SEEK_END);
    unsigned int lFileSize = ftell(pWavInfo->pHandle);
    fseek(pWavInfo->pHandle, 0, SEEK_SET);

    char rgchChunkID[4] = "\0\0\0";
    fread((void*)rgchChunkID, 4, 1, pWavInfo->pHandle);

    unsigned int lChunkSize = 0;
    fread((void*)&lChunkSize, 1, 4, pWavInfo->pHandle);

    char rgchFormat[4] = "\0\0\0";
    fread((void*)rgchFormat, 4, 1, pWavInfo->pHandle);

    char rgchSubchunk1ID[4] = "\0\0\0";
    fread((void*)rgchSubchunk1ID, 4, 1, pWavInfo->pHandle);

    unsigned int lSubchunk1Size = 0;
    fread((void*)&lSubchunk1Size, 1, 4, pWavInfo->pHandle);

    unsigned short wAudioFormat = 0;
    fread((void*)&wAudioFormat, 1, 2, pWavInfo->pHandle);

    fread((void*)&pWavInfo->wNumChannels, 1, 2, pWavInfo->pHandle);

    fread((void*)&pWavInfo->lSampleRate, 1, 4, pWavInfo->pHandle);

    unsigned int lByteRate = 0;
    fread((void*)&lByteRate, 1, 4, pWavInfo->pHandle);

    unsigned short wBlockAlign = 0;
    fread((void*)&wBlockAlign, 1, 2, pWavInfo->pHandle);

    unsigned short wBitsPerSample;
    fread((void*)&wBitsPerSample, 1, 2, pWavInfo->pHandle);

    char rgchSubchunk2ID[4] = "\0\0\0";
    fread((void*)rgchSubchunk2ID, 4, 1, pWavInfo->pHandle);
    
    unsigned int lSubchunk2Size = 0;
    fread((void*)&lSubchunk2Size, 1, 4, pWavInfo->pHandle);

    int result = 0;
    if (wBitsPerSample % 8 != 0)
        result = -1;

    pWavInfo->wBytesPerSample = wBitsPerSample / 8;
    if (pWavInfo->wBytesPerSample != 1 && pWavInfo->wBytesPerSample != 2 && pWavInfo->wBytesPerSample != 4)
        result = -2;

    unsigned int lFrameSize = pWavInfo->wNumChannels * pWavInfo->wBytesPerSample;
    if (lSubchunk2Size % lFrameSize != 0)
        result = -3;

    pWavInfo->lNumFrames = lSubchunk2Size / lFrameSize;
    pWavInfo->lNumDataBytes = lSubchunk2Size;

    return 0;
}

FILE* WaveFile::openWrite(char* pszPath)
{
    if (NULL == pszPath)
        return NULL;

    FILE* pHandle = fopen(pszPath, "wb");
    return pHandle;
}

int WaveFile::read(WavInfo* pWavInfo, float* pData)
{
    if (NULL == pWavInfo || NULL == pData)
        return -1;

    fseek(pWavInfo->pHandle, 44, SEEK_SET);
    unsigned char* pTemp = new unsigned char[pWavInfo->lNumDataBytes];
    fread(pTemp, 1, pWavInfo->lNumDataBytes, pWavInfo->pHandle);

    int lNumSamples = pWavInfo->lNumDataBytes / pWavInfo->wBytesPerSample;
    switch (pWavInfo->wBytesPerSample)
    {
        case 1:
        for (int i = 0; i < lNumSamples; i++)
            pData[i] = pTemp[i] / 127.0f - 1.0f;
        break;

        case 2:
        short* pShort = (short*)pTemp;
        for (int i = 0; i < lNumSamples; i++)
            pData[i] = pShort[i] / 32767.0f;
        break;

        case 4:
        int* pInt = (int*)pTemp;
        for (int i = 0; i < lNumSamples; i++)
            pData[i] = pInt[i] / 2147483647.0f;
        break;

        default:
        return -1;
    }

    delete[] pTemp;

    return 0;
}

int WaveFile::write(WavInfo* pWavInfo, float* pData)
{
    if (NULL == pWavInfo)
        return -1;

    if (NULL == pWavInfo->pHandle)
        return -2;

    union MultiPtr {
        void* pVoid;
        unsigned char* pUChar;
        short* pShort;
        int* pInt;
    };
    
    MultiPtr pTemp;
    pTemp.pVoid = NULL;
    int lNumSamples = pWavInfo->lNumDataBytes / pWavInfo->wBytesPerSample;
    switch (pWavInfo->wBytesPerSample)
    {
        case 1:
        pTemp.pUChar = new unsigned char[lNumSamples];
        for (int i = 0; i < lNumSamples; i++)
            pTemp.pUChar[i] = (unsigned char)((pData[i] + 1.0f) * 127.0f + 0.5f);
        break;

        case 2:
        pTemp.pShort = new short[lNumSamples];
        for (int i = 0; i < lNumSamples; i++)
            pTemp.pShort[i] = (short)((pData[i] + 1.0f) * 32767.0f + 0.5f) - 32767;
        break;

        case 4:
        pTemp.pInt = (void*)(new int[lNumSamples]);
        for (int i = 0; i < lNumSamples; i++)
            pTemp.pInt[i] = (int)((pData[i] + 1.0f) * 2147483647.0f + 0.5f) - 2147483647;
        break;

        default:
        return -3;
    }

    int itemCount = 0;
    
    itemCount += fwrite("RIFF", 1, 4, pWavInfo->pHandle);

    unsigned int lChunkSize = pWavInfo->lNumDataBytes - 8;
    itemCount += fwrite(&lChunkSize, 1, 4, pWavInfo->pHandle);

    itemCount += fwrite("WAVE", 1, 4, pWavInfo->pHandle);

    itemCount += fwrite("fmt ", 1, 4, pWavInfo->pHandle);

    unsigned int lSubchunk1Size = 16;
    itemCount += fwrite(&lSubchunk1Size, 1, 4, pWavInfo->pHandle);

    unsigned short wAudioFormat = 1;
    itemCount += fwrite(&wAudioFormat, 1, 2, pWavInfo->pHandle);

    itemCount += fwrite(&pWavInfo->wNumChannels, 1, 2, pWavInfo->pHandle);

    itemCount += fwrite(&pWavInfo->lSampleRate, 1, 4, pWavInfo->pHandle);

    unsigned short wBlockAlign = pWavInfo->wNumChannels * pWavInfo->wBytesPerSample;
    unsigned int lByteRate = pWavInfo->lSampleRate * wBlockAlign;
    itemCount += fwrite(&lByteRate, 1, 4, pWavInfo->pHandle);
    itemCount += fwrite(&wBlockAlign, 1, 2, pWavInfo->pHandle);

    unsigned short wBitsPerSample = pWavInfo->wBytesPerSample * 8;
    itemCount += fwrite(&wBitsPerSample, 1, 2, pWavInfo->pHandle);

    itemCount += fwrite("data", 1, 4, pWavInfo->pHandle);

    itemCount += fwrite(&pWavInfo->lNumDataBytes, 1, 4, pWavInfo->pHandle);

    switch (pWavInfo->wBytesPerSample)
    {
        case 1:
        itemCount += fwrite(pTemp.pUChar, 1, lNumSamples, pWavInfo->pHandle);
        delete[] pTemp.pUChar;
        break;

        case 2:
        itemCount += fwrite(pTemp.pShort, 2, lNumSamples, pWavInfo->pHandle);
        delete[] pTemp.pShort;
        break;

        case 4:
        itemCount += fwrite(pTemp.pInt, 4, lNumSamples, pWavInfo->pHandle);
        delete[] pTemp.pInt;
        break;

        default:
        return -4;
    }
    
    return itemCount;
}

int WaveFile::test()
{
    WavInfo wavInfo;
    
    gWaveFile.openRead("c:\\windows\\media\\chimes.wav", &wavInfo);
    float* pData = new float[wavInfo.lNumFrames * wavInfo.wNumChannels];
    gWaveFile.read(&wavInfo, pData);
    fclose(wavInfo.pHandle);

    wavInfo.pHandle = gWaveFile.openWrite("c:\\users\\owner\\desktop\\chimes2.wav");
    printf("write: %d\n", gWaveFile.write(&wavInfo, pData));
    wavInfo.pHandle = (FILE*)NULL;
    
    delete[] pData;
}


////////////////////////////////////////////////////////////////////////////////


// -1: Bad pWinData value
// -2: Internal allocation error
int asymWindow(double** pWinData, int nWinType, int nSectLen,
               int nNumAscSects, int nNumDescSects,
               int nNumZeroPreSects=0, int nNumZeroPostSects=0)
{
    // T is for type. It isn't used yet.
  
    if (NULL == pWinData)
        return -1;
    
    int nReturnLen = 0;
    
    int nWinLengthAsc  = nSectLen * nNumAscSects;
    int nWinLengthDesc = nSectLen * nNumDescSects;

    int nNumPreZeros  = nSectLen * nNumZeroPreSects;
    int nNumPostZeros = nSectLen * nNumZeroPostSects;
    
    nReturnLen = nNumPreZeros + nWinLengthAsc + nWinLengthDesc + nNumZeroPostSects;
    *pWinData = new double[nReturnLen];
    if (NULL == *pWinData)
        return -2;

    double* pConveniencePtr = *pWinData;
    memset(pConveniencePtr, 0, nReturnLen);
    int nAscSectStart  = nNumPreZeros;
    int nDescSectStart = nAscSectStart  + nWinLengthAsc;
    int nDescSectEnd   = nDescSectStart + nWinLengthDesc;
    
    for (int i = nAscSectStart; i < nDescSectStart; i++)
    {
        int nAscIndex = i - nAscSectStart;
        pConveniencePtr[i] = (double)nAscIndex / (double)nWinLengthAsc;
    }
    
    for (int i = nDescSectStart; i < nDescSectEnd; i++)
    {
        int nDescIndex = i - nDescSectStart;
        pConveniencePtr[i] = 1.0f - (double)nDescIndex / (double)nWinLengthDesc;
    }
    
    gMutagenTestGraph = new TGraph(nReturnLen, (const double*)gMutagenOrdinates, 
                                   (const double*)pConveniencePtr);
    gMutagenTestGraph->Draw("AL");
    
    return nReturnLen;
}
/*
void asymWindowTest()
{
    if (NULL != gMutagenTestGraph) {
        delete gMutagenTestGraph;
        gMutagenTestGraph = NULL;
    }
    
    float* pWinData;
    int result = asymWindow(&pWinData, 0, 256, 1, 3);
    print("result = $d\n", result);

    if (result > 0) {
        gMutagenTestGraph = new TGraph
    }
}
*/

//////////////////////////////////////////////////////////////////////////////

/*****************************************************************************
*                                                                            *
*       DIGITAL SIGNAL PROCESSING TOOLS                                      *
*       Version 1.01, 1999/11/07                                             *
*       (c) 1999 - Laurent de Soras                                          *
*                                                                            *
*       FFTReal.h                                                            *
*       Fourier transformation of real number arrays.                        *
*       Portable ISO C++                                                     *
*                                                                            *
* Tab = 3                                                                    *
*****************************************************************************/

class FFTReal
{

public:

    typedef double flt_t;

    explicit            FFTReal (const long length);
                        ~FFTReal ();
    void                do_fft (flt_t f [], const flt_t x []) const;
    void                do_ifft (const flt_t f [], flt_t x []) const;
    void                rescale (flt_t x []) const;


private:

    /* Bit-reversed look-up table nested class */
    class BitReversedLUT
    {
    public:
        explicit            BitReversedLUT (const int nbr_bits);
                            ~BitReversedLUT ();
        const long *    get_ptr () const
        {
            return (_ptr);
        }
    private:
        long *            _ptr;
    };

    /* Trigonometric look-up table nested class */
    class    TrigoLUT
    {
    public:
        explicit            TrigoLUT (const int nbr_bits);
                            ~TrigoLUT ();
        const flt_t    *    get_ptr (const int level) const
        {
            return (_ptr + (1L << (level - 1)) - 4);
        };
    private:
        flt_t    *            _ptr;
    };

    const BitReversedLUT    _bit_rev_lut;
    const TrigoLUT    _trigo_lut;
    const flt_t        _sqrt2_2;
    const long        _length;
    const int        _nbr_bits;
    flt_t *            _buffer_ptr;



/*\\\ FORBIDDEN MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

                       FFTReal (const FFTReal &other);
    const FFTReal&     operator = (const FFTReal &other);
    int                operator == (const FFTReal &other);
    int                operator != (const FFTReal &other);
};


/*****************************************************************************

    LEGAL

    Source code may be freely used for any purpose, including commercial
    applications. Programs must display in their "About" dialog-box (or
    documentation) a text telling they use these routines by Laurent de Soras.
    Modified source code can be distributed, but modifications must be clearly
    indicated.

    CONTACT

    Laurent de Soras
    92 avenue Albert 1er
    92500 Rueil-Malmaison
    France

    ldesoras@club-internet.fr

*****************************************************************************/

/*****************************************************************************
*                                                                            *
*       DIGITAL SIGNAL PROCESSING TOOLS                                      *
*       Version 1.01, 1999/11/07                                             *
*       (c) 1999 - Laurent de Soras                                          *
*                                                                            *
*       FFTReal.cpp                                                          *
*       Fourier transformation of real number arrays.                        *
*       Portable ISO C++                                                     *
*                                                                            *
* Tab = 3                                                                    *
*****************************************************************************/



/*\\\ INCLUDE FILES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/




/*\\\ PUBLIC MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/



/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - length: length of the array on which we want to do a FFT.       */
/*                  Range: power of 2 only, > 0.                            */
/*      Throws: std::bad_alloc, anything                                    */
/*==========================================================================*/

FFTReal::FFTReal (const long length)
:    _length (length)
,    _nbr_bits (int (floor (log ((double)length) / log (2.0) + 0.5)))
,    _bit_rev_lut (int (floor (log ((double)length) / log (2.0) + 0.5)))
,    _trigo_lut (int (floor (log ((double)length) / log (2.0) + 0.5)))
,    _sqrt2_2 (flt_t (sqrt (2.0) * 0.5))
{
    assert ((1L << _nbr_bits) == length);

    _buffer_ptr = 0;
    if (_nbr_bits > 2)
    {
        _buffer_ptr = new flt_t [_length];
    }
}



/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/

FFTReal::~FFTReal (void)
{
    delete [] _buffer_ptr;
    _buffer_ptr = 0;
}



/*==========================================================================*/
/*      Name: do_fft                                                        */
/*      Description: Compute the FFT of the array.                          */
/*      Input parameters:                                                   */
/*        - x: pointer on the source array (time).                          */
/*      Output parameters:                                                  */
/*        - f: pointer on the destination array (frequencies).              */
/*             f [0...length(x)/2] = real values,                           */
/*             f [length(x)/2+1...length(x)-1] = imaginary values of        */
/*               coefficents 1...length(x)/2-1.                             */
/*      Throws: Nothing                                                     */
/*==========================================================================*/

void    FFTReal::do_fft (flt_t f [], const flt_t x []) const
{

/*______________________________________________
 *
 * General case
 *______________________________________________
 */

    if (_nbr_bits > 2)
    {
        flt_t *            sf;
        flt_t *            df;

        if (_nbr_bits & 1)
        {
            df = _buffer_ptr;
            sf = f;
        }
        else
        {
            df = f;
            sf = _buffer_ptr;
        }

        /* Do the transformation in several pass */
        {
            int        pass;
            long        nbr_coef;
            long        h_nbr_coef;
            long        d_nbr_coef;
            long        coef_index;

            /* First and second pass at once */
            {
                const long* const    bit_rev_lut_ptr = _bit_rev_lut.get_ptr ();
                coef_index = 0;
                do
                {
                    const long        rev_index_0 = bit_rev_lut_ptr [coef_index];
                    const long        rev_index_1 = bit_rev_lut_ptr [coef_index + 1];
                    const long        rev_index_2 = bit_rev_lut_ptr [coef_index + 2];
                    const long        rev_index_3 = bit_rev_lut_ptr [coef_index + 3];

                    flt_t    * const    df2 = df + coef_index;
                    df2 [1] = x [rev_index_0] - x [rev_index_1];
                    df2 [3] = x [rev_index_2] - x [rev_index_3];

                    const flt_t        sf_0 = x [rev_index_0] + x [rev_index_1];
                    const flt_t        sf_2 = x [rev_index_2] + x [rev_index_3];

                    df2 [0] = sf_0 + sf_2;
                    df2 [2] = sf_0 - sf_2;
                    
                    coef_index += 4;
                }
                while (coef_index < _length);
            }

            /* Third pass */
            {
                coef_index = 0;
                const flt_t        sqrt2_2 = _sqrt2_2;
                do
                {
                    flt_t                v;

                    sf [coef_index] = df [coef_index] + df [coef_index + 4];
                    sf [coef_index + 4] = df [coef_index] - df [coef_index + 4];
                    sf [coef_index + 2] = df [coef_index + 2];
                    sf [coef_index + 6] = df [coef_index + 6];

                    v = (df [coef_index + 5] - df [coef_index + 7]) * sqrt2_2;
                    sf [coef_index + 1] = df [coef_index + 1] + v;
                    sf [coef_index + 3] = df [coef_index + 1] - v;

                    v = (df [coef_index + 5] + df [coef_index + 7]) * sqrt2_2;
                    sf [coef_index + 5] = v + df [coef_index + 3];
                    sf [coef_index + 7] = v - df [coef_index + 3];

                    coef_index += 8;
                }
                while (coef_index < _length);
            }

            /* Next pass */
            for (pass = 3; pass < _nbr_bits; ++pass)
            {
                coef_index = 0;
                nbr_coef = 1 << pass;
                h_nbr_coef = nbr_coef >> 1;
                d_nbr_coef = nbr_coef << 1;
                const flt_t    * const    cos_ptr = _trigo_lut.get_ptr (pass);
                do
                {
                    long                i;
                    const flt_t    *    const sf1r = sf + coef_index;
                    const flt_t    *    const sf2r = sf1r + nbr_coef;
                    flt_t *            const dfr = df + coef_index;
                    flt_t *            const dfi = dfr + nbr_coef;

                    /* Extreme coefficients are always real */
                    dfr [0] = sf1r [0] + sf2r [0];
                    dfi [0] = sf1r [0] - sf2r [0];    // dfr [nbr_coef] =
                    dfr [h_nbr_coef] = sf1r [h_nbr_coef];
                    dfi [h_nbr_coef] = sf2r [h_nbr_coef];

                    /* Others are conjugate complex numbers */
                    const flt_t    * const    sf1i = sf1r + h_nbr_coef;
                    const flt_t    * const    sf2i = sf1i + nbr_coef;
                    for (i = 1; i < h_nbr_coef; ++ i)
                    {
                        const flt_t        c = cos_ptr [i];                    // cos (i*PI/nbr_coef);
                        const flt_t        s = cos_ptr [h_nbr_coef - i];    // sin (i*PI/nbr_coef);
                        flt_t                v;

                        v = sf2r [i] * c - sf2i [i] * s;
                        dfr [i] = sf1r [i] + v;
                        dfi [-i] = sf1r [i] - v;    // dfr [nbr_coef - i] =

                        v = sf2r [i] * s + sf2i [i] * c;
                        dfi [i] = v + sf1i [i];
                        dfi [nbr_coef - i] = v - sf1i [i];
                    }

                    coef_index += d_nbr_coef;
                }
                while (coef_index < _length);

                /* Prepare to the next pass */
                {
                    flt_t    * const        temp_ptr = df;
                    df = sf;
                    sf = temp_ptr;
                }
            }
        }
    }

/*______________________________________________
 *
 * Special cases
 *______________________________________________
 */

    /* 4-point FFT */
    else if (_nbr_bits == 2)
    {
        f [1] = x [0] - x [2];
        f [3] = x [1] - x [3];

        const flt_t            b_0 = x [0] + x [2];
        const flt_t            b_2 = x [1] + x [3];
        
        f [0] = b_0 + b_2;
        f [2] = b_0 - b_2;
    }

    /* 2-point FFT */
    else if (_nbr_bits == 1)
    {
        f [0] = x [0] + x [1];
        f [1] = x [0] - x [1];
    }

    /* 1-point FFT */
    else
    {
        f [0] = x [0];
    }
}



/*==========================================================================*/
/*      Name: do_ifft                                                       */
/*      Description: Compute the inverse FFT of the array. Notice that      */
/*                   IFFT (FFT (x)) = x * length (x). Data must be          */
/*                   post-scaled.                                           */
/*      Input parameters:                                                   */
/*        - f: pointer on the source array (frequencies).                   */
/*             f [0...length(x)/2] = real values,                           */
/*             f [length(x)/2+1...length(x)] = imaginary values of          */
/*               coefficents 1...length(x)-1.                               */
/*      Output parameters:                                                  */
/*        - x: pointer on the destination array (time).                     */
/*      Throws: Nothing                                                     */
/*==========================================================================*/

void    FFTReal::do_ifft (const flt_t f [], flt_t x []) const
{

/*______________________________________________
 *
 * General case
 *______________________________________________
 */

    if (_nbr_bits > 2)
    {
        flt_t *            sf = const_cast <flt_t *> (f);
        flt_t *            df;
        flt_t *            df_temp;

        if (_nbr_bits & 1)
        {
            df = _buffer_ptr;
            df_temp = x;
        }
        else
        {
            df = x;
            df_temp = _buffer_ptr;
        }

        /* Do the transformation in several pass */
        {
            int            pass;
            long            nbr_coef;
            long            h_nbr_coef;
            long            d_nbr_coef;
            long            coef_index;

            /* First pass */
            for (pass = _nbr_bits - 1; pass >= 3; --pass)
            {
                coef_index = 0;
                nbr_coef = 1 << pass;
                h_nbr_coef = nbr_coef >> 1;
                d_nbr_coef = nbr_coef << 1;
                const flt_t    *const cos_ptr = _trigo_lut.get_ptr (pass);
                do
                {
                    long                i;
                    const flt_t    *    const sfr = sf + coef_index;
                    const flt_t    *    const sfi = sfr + nbr_coef;
                    flt_t *            const df1r = df + coef_index;
                    flt_t *            const df2r = df1r + nbr_coef;

                    /* Extreme coefficients are always real */
                    df1r [0] = sfr [0] + sfi [0];        // + sfr [nbr_coef]
                    df2r [0] = sfr [0] - sfi [0];        // - sfr [nbr_coef]
                    df1r [h_nbr_coef] = sfr [h_nbr_coef] * 2;
                    df2r [h_nbr_coef] = sfi [h_nbr_coef] * 2;

                    /* Others are conjugate complex numbers */
                    flt_t * const    df1i = df1r + h_nbr_coef;
                    flt_t * const    df2i = df1i + nbr_coef;
                    for (i = 1; i < h_nbr_coef; ++ i)
                    {
                        df1r [i] = sfr [i] + sfi [-i];        // + sfr [nbr_coef - i]
                        df1i [i] = sfi [i] - sfi [nbr_coef - i];

                        const flt_t        c = cos_ptr [i];                    // cos (i*PI/nbr_coef);
                        const flt_t        s = cos_ptr [h_nbr_coef - i];    // sin (i*PI/nbr_coef);
                        const flt_t        vr = sfr [i] - sfi [-i];        // - sfr [nbr_coef - i]
                        const flt_t        vi = sfi [i] + sfi [nbr_coef - i];

                        df2r [i] = vr * c + vi * s;
                        df2i [i] = vi * c - vr * s;
                    }

                    coef_index += d_nbr_coef;
                }
                while (coef_index < _length);

                /* Prepare to the next pass */
                if (pass < _nbr_bits - 1)
                {
                    flt_t    * const    temp_ptr = df;
                    df = sf;
                    sf = temp_ptr;
                }
                else
                {
                    sf = df;
                    df = df_temp;
                }
            }

            /* Antepenultimate pass */
            {
                const flt_t        sqrt2_2 = _sqrt2_2;
                coef_index = 0;
                do
                {
                    df [coef_index] = sf [coef_index] + sf [coef_index + 4];
                    df [coef_index + 4] = sf [coef_index] - sf [coef_index + 4];
                    df [coef_index + 2] = sf [coef_index + 2] * 2;
                    df [coef_index + 6] = sf [coef_index + 6] * 2;

                    df [coef_index + 1] = sf [coef_index + 1] + sf [coef_index + 3];
                    df [coef_index + 3] = sf [coef_index + 5] - sf [coef_index + 7];

                    const flt_t        vr = sf [coef_index + 1] - sf [coef_index + 3];
                    const flt_t        vi = sf [coef_index + 5] + sf [coef_index + 7];

                    df [coef_index + 5] = (vr + vi) * sqrt2_2;
                    df [coef_index + 7] = (vi - vr) * sqrt2_2;

                    coef_index += 8;
                }
                while (coef_index < _length);
            }

            /* Penultimate and last pass at once */
            {
                coef_index = 0;
                const long *    bit_rev_lut_ptr = _bit_rev_lut.get_ptr ();
                const flt_t    *    sf2 = df;
                do
                {
                    {
                        const flt_t        b_0 = sf2 [0] + sf2 [2];
                        const flt_t        b_2 = sf2 [0] - sf2 [2];
                        const flt_t        b_1 = sf2 [1] * 2;
                        const flt_t        b_3 = sf2 [3] * 2;

                        x [bit_rev_lut_ptr [0]] = b_0 + b_1;
                        x [bit_rev_lut_ptr [1]] = b_0 - b_1;
                        x [bit_rev_lut_ptr [2]] = b_2 + b_3;
                        x [bit_rev_lut_ptr [3]] = b_2 - b_3;
                    }
                    {
                        const flt_t        b_0 = sf2 [4] + sf2 [6];
                        const flt_t        b_2 = sf2 [4] - sf2 [6];
                        const flt_t        b_1 = sf2 [5] * 2;
                        const flt_t        b_3 = sf2 [7] * 2;

                        x [bit_rev_lut_ptr [4]] = b_0 + b_1;
                        x [bit_rev_lut_ptr [5]] = b_0 - b_1;
                        x [bit_rev_lut_ptr [6]] = b_2 + b_3;
                        x [bit_rev_lut_ptr [7]] = b_2 - b_3;
                    }

                    sf2 += 8;
                    coef_index += 8;
                    bit_rev_lut_ptr += 8;
                }
                while (coef_index < _length);
            }
        }
    }

/*______________________________________________
 *
 * Special cases
 *______________________________________________
 */

    /* 4-point IFFT */
    else if (_nbr_bits == 2)
    {
        const flt_t        b_0 = f [0] + f [2];
        const flt_t        b_2 = f [0] - f [2];

        x [0] = b_0 + f [1] * 2;
        x [2] = b_0 - f [1] * 2;
        x [1] = b_2 + f [3] * 2;
        x [3] = b_2 - f [3] * 2;
    }

    /* 2-point IFFT */
    else if (_nbr_bits == 1)
    {
        x [0] = f [0] + f [1];
        x [1] = f [0] - f [1];
    }

    /* 1-point IFFT */
    else
    {
        x [0] = f [0];
    }
}



/*==========================================================================*/
/*      Name: rescale                                                       */
/*      Description: Scale an array by divide each element by its length.   */
/*                   This function should be called after FFT + IFFT.       */
/*      Input/Output parameters:                                            */
/*        - x: pointer on array to rescale (time or frequency).             */
/*      Throws: Nothing                                                     */
/*==========================================================================*/

void    FFTReal::rescale (flt_t x []) const
{
    const flt_t        mul = flt_t (1.0 / _length);
    long                i = _length - 1;

    do
    {
        x [i] *= mul;
        --i;
    }
    while (i >= 0);
}



/*\\\ NESTED CLASS MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/



/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - nbr_bits: number of bits of the array on which we want to do a  */
/*                    FFT. Range: > 0                                       */
/*      Throws: std::bad_alloc                                              */
/*==========================================================================*/

FFTReal::BitReversedLUT::BitReversedLUT (const int nbr_bits)
{
    long                length;
    long                cnt;
    long                br_index;
    long                bit;

    length = 1L << nbr_bits;
    _ptr = new long [length];

    br_index = 0;
    _ptr [0] = 0;
    for (cnt = 1; cnt < length; ++cnt)
    {
        /* ++br_index (bit reversed) */
        bit = length >> 1;
        while (((br_index ^= bit) & bit) == 0)
        {
            bit >>= 1;
        }

        _ptr [cnt] = br_index;
    }
}



/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/

FFTReal::BitReversedLUT::~BitReversedLUT (void)
{
    delete [] _ptr;
    _ptr = 0;
}



/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - nbr_bits: number of bits of the array on which we want to do a  */
/*                    FFT. Range: > 0                                       */
/*      Throws: std::bad_alloc, anything                                    */
/*==========================================================================*/

FFTReal::TrigoLUT::TrigoLUT (const int nbr_bits)
{
    long        total_len;

    _ptr = 0;
    if (nbr_bits > 3)
    {
        total_len = (1L << (nbr_bits - 1)) - 4;
        _ptr = new flt_t [total_len];

        const double    PI = atan (1.0) * 4;
        for (int level = 3; level < nbr_bits; ++level)
        {
            const long        level_len = 1L << (level - 1);
            flt_t    * const    level_ptr = const_cast<flt_t *> (get_ptr (level));
            const double    mul = PI / (level_len << 1);

            for (long i = 0; i < level_len; ++ i)
            {
                level_ptr [i] = (flt_t) cos (i * mul);
            }
        }
    }
}



/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/

FFTReal::TrigoLUT::~TrigoLUT (void)
{
    delete [] _ptr;
    _ptr = 0;
}




/*****************************************************************************

    LEGAL

    Source code may be freely used for any purpose, including commercial
    applications. Programs must display in their "About" dialog-box (or
    documentation) a text telling they use these routines by Laurent de Soras.
    Modified source code can be distributed, but modifications must be clearly
    indicated.

    CONTACT

    Laurent de Soras
    92 avenue Albert 1er
    92500 Rueil-Malmaison
    France

    ldesoras@club-internet.fr

*****************************************************************************/



class STFT {
public:
    STFT(int nInputLen, double* pfInputReal, double* pfInputCplx, int nWinLen,
         double* pfWindow, int nOverlapNumber)
    :    fTransform(nWinLen), fInputReal(pfInputReal),  fInputCplx(pfInputCplx),
         fInputLen(nInputLen), fWindow(pfWindow), fOverlapNumber(nOverlapNumber),
         fNumSpectra(0), fSpectraReal(NULL), fSpectraCplx(NULL),
         fSpectraRealFlattened(NULL), fSpectraCplxFlattened(NULL) {}
    int Run();
    
    int fOverlapNumber;
    int fNumSpectra;
    double** fSpectraReal;
    double** fSpectraCplx;
    
private:
    FFTReal* fTransform;
    double* fInputReal;
    double* fInputCplx;
    int fInputLen;
    double* fWindow;
    double* fSpectraRealFlattened;
    double* fSpectraCplxFlattened;
};

int STFT::Run()
{
    if (0 == fOverlapNumber)
        return -1;
        
    int nStrideLen = fWinLen / fOverlapNumber;
    // If fOverlapNumber doesn't divide fWinLen
    if (nStrideLen * fOverlapNumber != fWinLen)
        return -2;

    if (NULL == fTransform)
        return -3;
        
    if (NULL == fInputReal || NULL == fInputCplx)
        return -4;
        
    if (0 >= fInputLen)
        return 0;
        
    int numStrides = (fInputLen-1) / nStrideLen + 1;
    if (0 >= numStrides)
        return 0;
        
    double* pfFFTBufferReal = new double[nStrideLen];
    double* pfFFTBufferCplx = new double[nStrideLen];
    fSpectraRealFlattened = new double[fTransform->GetSize() * numStrides];
    fSpectraCplxFlattened = new double[fTransform->GetSize() * numStrides];
    fSpectraReal = new (double*)[numStrides];
    fSpectraCplx = new (double*)[numStrides];
    if (NULL == pfFFTBufferReal || NULL == pfFFTBufferCplx ||
        NULL == fSpectraRealFlattened || NULL == fSpectraCplxFlattened ||
        NULL == fSpectraReal || NULL == fSpectraCplx)
    {
        mutagenSafeDelete(pfFFTBufferReal);
        mutagenSafeDelete(pfFFTBufferCplx);
        mutagenSafeDelete(fSpectraRealFlattened);
        mutagenSafeDelete(fSpectraCplxFlattened);
        mutagenSafeDelete(fSpectraReal);
        mutagenSafeDelete(fSpectraCplx);
        return -5;
    }
    
    for (int i = 0; i < numStrides; i++) {
        fSpectraReal[i] = fSpectraRealFlattened + i*nStrideLen;
        fSpectraCplx[i] = fSpectraCplxFlattened + i*nStrideLen;
    }
    
    int nMaxEvenLength = (fInputLen / nStrideLen) * nStrideLen;    
    for (int i = 0; i < nMaxEvenLength; i += nStrideLen) {
        int nStrideNumber = i / nStrideLen;
        memcpy(pfFFTBufferReal, fInput+i, fTransform->GetSize());
        memset(pfFFTBufferCplx, 0, fTransform->GetSize());
        fTransform->SetPoints((const double*)pfFFTBufferReal,
                              (const double*)pfFFTBufferCplx);
        fTransform->Transform();
        fTransform->GetPoints(pfFFTBufferReal, pfFFTBufferCplx);
        memcpy(fSpectraReal[i], pfFFTBufferReal, fTransform->GetSize());
        memcpy(fSpectraCplx[i], pfFFTBufferCplx, fTransform->GetSize());
    }
    
    mutagenSafeDelete(pfFFTBufferReal, true);
    mutagenSafeDelete(pfFFTBufferCplx, true);
    mutagenSafeDelete(fSpectraRealFlattened, true);
    mutagenSafeDelete(fSpectraCplxFlattened, true);
    mutagenSafeDelete(fSpectraReal, true);
    mutagenSafeDelete(fSpectraCplx, true);
    
    fNumSpectra = numStrides;
    return numStrides;
}
/*
function out = stftbrent(x, w, overlapnumber)

  x = x(:);
  w = w(:);
  
  xlen = size(x, 1);
  winlength = size(w, 1);
  stride = winlength / overlapnumber;
  numstrides = (floor(xlen/stride) + overlapnumber - 1);
  fullxlen = numstrides * stride;
  fullx = [x ; zeros(fullxlen-xlen, 1)];

  y = zeros(winlength, numstrides);
  
  for winstart = 1:stride:(fullxlen-winlength)  
    n = (winstart - 1) / stride + 1;  
    winend = winstart + winlength - 1;  
    fftin = w .* fullx(winstart:winend);  
    y(:, n) = fft(fftin);
  endfor  
  
  out = y / winlength;  

endfunction
*/

TCanvas *c1;

void rose_image()
{
   // Display image in a new canvas and pad.
   //Author: Valeriy Onuchin
   
   TImage *img = TImage::Open("image/rose512.jpg");

   if (!img) {
      printf("Could not create an image... exit\n");
      return;
   }

   img->SetConstRatio(0);
   img->SetImageQuality(TAttImage::kImgBest);

   TString fp = gEnv->GetValue("Root.TTFontPath", "");
   TString bc = fp + "/BlackChancery.ttf";
   TString ar = fp + "/arial.ttf";

   // draw text over image with funny font
   img->DrawText(120, 160, "Hello World!", 32, 
                 gROOT->GetColor(4)->AsHexString(), 
                 bc, TImage::kShadeBelow);

   // draw text over image with foreground specified by pixmap
   img->DrawText(250, 350, "goodbye cruel world ...", 24, 0, 
                 ar, TImage::kPlain, "fore.xpm");

   TImage *img2 = TImage::Open("image/mditestbg.xpm");

   // tile image
   img2->Tile(img->GetWidth(), img->GetHeight());

   c1 = new TCanvas("poop", "examples of image manipulations", 760, 900);
   c1->Divide(2, 3);
   c1->cd(1);
   img->Draw("xxx");
   img->SetEditable(kTRUE);

   c1->cd(2);
   // averaging with mditestbg.xpm image
   TImage *img3 = (TImage*)img->Clone("img3");
   img3->Merge(img2, "allanon");
   img3->Draw();

   // contrasting (tint with itself)
   c1->cd(3);
   TImage *img4 = (TImage*)img->Clone("img4");
   img4->Merge(img4, "tint");

   // draw filled rectangle with magenta color
   img4->FillRectangle("#FF00FF", 20, 220, 40, 40);

   // Render multipoint alpha-blended gradient (R->G->B)
   img4->Gradient(0, "#FF0000 #00FF00 #220000FF", 0, 50, 50, 100, 100);

   // draw semi-transparent 3D button
   img4->Bevel(300, 20, 160, 40, "#ffffffff", "#fe000000", 3, 0);
   img4->DrawLine(10, 100, 100, 10, "#0000ff", 4);
   img4->Draw();

   // vectorize image. Reduce palette to 256 colors
   c1->cd(4);
   TImage *img5 = (TImage*)img->Clone("img5");
   img5->Vectorize(256);
   img5->Draw();

   // quantization of the image
   c1->cd(5);
   TImage *img6 = (TImage*)img->Clone("img6");
   TImagePalette *pal = (TImagePalette *)&img5->GetPalette();
   TArrayD *arr = img6->GetArray(50, 40, pal);
   img6->SetImage(arr->GetArray(), 50, 40, pal);
   img6->Draw();

   // HSV adjustment (convert red to yellow)
   c1->cd(6);
   TImage *img7 = (TImage*)img->Clone("img7");
   img7->HSV(0, 40, 40);
   img7->Draw();

}