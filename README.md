# Design_of_DNA_based_Color_Images_Cryptosystem
## Introduction
DNA based color images refer to the encoding of the color images into DNA sequences called as DNA computing. 
DNA computing is a cutting-edge technology that harnesses the power of genetic material to perform complex computations. 
By encoding information into DNA strands, scientists can create vast amounts of data storage in a space smaller than a single grain of sand. 
The basic idea behind DNA computing is that DNA can be used to store and process information in a way that is similar to how a computer uses binary digits (bits) to represent data. 
Instead of representing images using traditional pixel values the image data is converted into a sequence of DNA nucleotides(Adenine-A, Thymine-T, Cytosine-C and Guanine-G). 
The process of converting color images’s pixel values into a DNA sequence involve mapping the pixel value of a image to specific sequence of DNA bases. 
For each color channel i.e. RGB of pixel might be represented by a certain DNA bases or combination of bases. 
Overview of Conversion- The pixel value 0-255 can be converted into DNA sequence. Using a set of predefined rules we can map  0-63 to A, 64-127 to T, 128-191 to C and 192-255 to G.   

## The Proposed Cryptosystem
The cryptosystem makes use of six 2D-Multiple Collapse Chaotic Maps, which are produced by the Key Scheduling Algorithm from a 240-symmetric key. The first stage involves intermixing of channels of the image then utilising the proposed operator ⊗ to operate the rows of the intermixed images with the initial array. Using a modified version of Arnold's Cat Map [10], which can be applied to both square and rectangular images, the pixels of the composite image are jumbled. The suggested operator is then used to further mix the pixels in the directions of forward spiral row, forward spiral column, reverse spiral row, and reverse spiral column. These keyless procedures are necessary to eliminate the link between pixels and the image's visual perception. Lastly, a non-linear layer has been added by combining the suggested operator with the idea of DNA cryptography. One of the eight encoding principle identified by the rule chart is applied to the mixed image that was obtained in the preceding phases, converting it into a sequence of DNA nucleotides. At this stage DNA encoded image is obtained and these encodings are further substituted using proposed operator with the substitution map. Lastly, one of the eight DNA decoding guidelines identified by a different rule chart is applied to decode the substituted DNA sequences. 

## Encryption Process
1. Inter Channel Mixing
2. Mix rows with Initial Vector
3. Arnold Cat Map algorithm
4. Spiral Mixing in four different directions
5. DNA Encoding
6. DNA Substitution
7. DNA Decoding
8. Final cipher image

## Decryption
Reverse of Encryption procedure.

## Experimental Results and Analysis
### Metrics
1. Key Space and Key Sensitivity.
2. Histogram Analysis
3. Variance Analysis
4. chi Square Analysis
5. Entropy Analysis
6. Correlation Analysis
7. NPCR and UACI
8. PSNR
9. Comparative Analysis with other operator and Referenced Present Cryptosystems.
