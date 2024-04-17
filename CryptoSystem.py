#!/usr/bin/env python
# coding: utf-8

# In[5]:


import math
import secrets
import hashlib
import hmac
from decimal import Decimal, getcontext

class Cryptosystem:
    def __init__(self, image):
        self.m = image.shape[0]
        self.n = image.shape[1]
        print('m,n:',self.m,self.n)
        
    def generate_symmetric_key(self, key_length_bits):
        random_bytes = secrets.token_bytes(key_length_bits // 8)  # Convert bits to bytes
        return random_bytes

    def derive_symmetric_key(self, secret, salt, info, key_length):
        hmac_object = hmac.new(salt, secret, hashlib.sha256)
        derived_key = hmac_object.digest()
        if len(derived_key) < key_length:
            derived_key += hmac_object.digest(secret) 
        return derived_key[:key_length]
    
    def getSymmetricKey(self, length):
        symmetric_key  = self.generate_symmetric_key(length)
        salt = b"SomeRandomSalt" 
        info = b"AdditionalInfo"
        final_symmetric_key = self.derive_symmetric_key(symmetric_key, salt, info, length // 8)  # Convert bits to bytes
        secure_key = final_symmetric_key.hex()
        print("Symmetric Key (hexadecimal):", secure_key, 'with length:', len(secure_key))
        binary_key = bin(int(secure_key, 16))[2:]
        print('binary key: \n',binary_key,'\n length of key:', len(binary_key))
        return binary_key
        
    
    def generateInitials(self, binary_key):
        a0=b0=x=y=t=c1=c2=0
        for i in range(0,40):
            a0 += (int(binary_key[i]) * 2**(40-i))/(2**40)

        for i in range(40,80):
            b0 += (int(binary_key[i]) * 2**(80-i))/(2**40)

        for i in range(80,120):
            x += (int(binary_key[i]) * 2**(120-i))/(2**40)

        for i in range(120,160):
            y += (int(binary_key[i]) * 2**(160-i))/(2**40)

        for i in range(160,200):
            t += (int(binary_key[i]) * 2**(200-i))/(2**40)

        for i in range(200,220):
            c1 += (int(binary_key[i]) * 2**(220-i))/(2**20)

        for i in range(220,240):
            c2 += (int(binary_key[i]) * 2**(240-i))/(2**20)
            
        a = (a0 + t * c1) % 5 + 16 
        b = (b0 + t * c2) % 5 + 16
        x0 = (x + t * c1) % 2 - 1
        y0 = (y + t * c2) % 2 - 1
        
        return a,b,x0,y0
    
    def twoDMCCM(self, a, b, x, y):
        x1 = math.atan((b / (10 * a * y)) + math.tan(a * math.pi * x))
        y1 = math.atan((b / (a * x)) + math.tan(10 *a * math.pi * y ))
        return x1,y1
    
    def getInitialVector(self, pd, a, b, x0, y0):
        IV = []
        for i in range(0, self.m):
            x1,y1 = self.twoDMCCM(a, b, x0, y0)
            IV.append(abs(x1) * abs(y1) * ((pow(2,31)-1) % pow(2,pd)))
            x0 = x1
            y0 = y1
        return IV
    
    def rm_freq(self, RM):
        freq={}
        for y in RM:
            for z in y:
                if z in freq:
                    freq[z] +=1
                else:
                    freq[z] = 1
        #print(freq)
        myKeys = list(freq.keys())
        myKeys.sort()
        sorted_dict = {i: freq[i] for i in myKeys}

        print("No of times rule present in rule map")
        print(sorted_dict)
    
    def getRuleMap(self,a,b,pd):
        newp= (int)(pd/4)
        RM = []
        for i in range(1,self.m + 1):
            count=1 
            for j in range(1,self.n + 1):
                pos=(i-1)*self.n+j
                Xpos,Ypos = self.twoDMCCM(a,b,pos,pos)
                Ux = abs(Xpos) * ((pow(2,31)-1) % pow(2,pd))
                Uy = abs(Ypos) * ((pow(2,31)-1) % pow(2,pd))
                lit=[]
                for k in range(0,newp):
                    lit.insert(count,int(Ux % 8 +1))
                    lit.insert(count+1,int(Uy % 8 + 1)) 
                    count = count + 2;
                    Ux = math.floor(Ux/8)
                    Uy = math.floor(Uy/8)
                RM.append(lit)
        self.rm_freq(RM)
        return RM
    
    EncodingRules = {
        "00" : ['A','A','C','G','C','G','T','T'],
        "01" : ['C','G','A','A','T','T','C','G'],
        "10" : ['G','C','T','T','A','A','G','C'],
        "11" : ['T','T','G','C','G','C','A','A'],
    }
    
    def binaryEquivalent(self, U, pd): 
        getcontext().prec = pd
        d = Decimal(U)
        i = int(d * (1 << pd-1))
        bn = bin(i).replace("0b", "")
        if(len(bn)>pd):
            return bn[0:pd]
        else:
            for i in range(1,pd-len(bn)+1):
                bn='0'+bn     
        return bn
    
    def DSM_freq(self, DSM):
        all_freq= {}
        for x in DSM:
            if x in all_freq:
                all_freq[x] += 1
            else:
                all_freq[x] = 1
        print("No of times Neucleotides present in DSMap")
        print(all_freq)
    
    def getDSMMap(self,a,b, RM):
        DSM = []
        pd=8
        for i in range(1,self.m+1):
            count = 0
            for j in range(1,self.n+1):
                pos = (i-1)*self.n+j
                Xpos,Ypos = self.twoDMCCM(a,b,pos,pos)
                U = abs(Xpos) * abs(Ypos) * ((pow(2,31)-1) % pow(2,pd))
                bVal = self.binaryEquivalent(U, pd) #binary equivalent of U in p bit
                lst=[]
                for k in range(0,len(bVal),2):
                    str=bVal[k]+bVal[k+1]
                    lit = self.EncodingRules.get(str)
                    DSM.append(lit[RM[((i-1)*self.n+(j-1))][count%4]-1])
                    count = count + 1
        self.DSM_freq(DSM)
        return DSM
    

    


# In[12]:


import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
class Encryption(Cryptosystem):
    def __init__(self, image):
        self.m = image.shape[0]
        self.n = image.shape[1]
        blue, green, red = cv2.split(image)
        self.red = red
        self.green = green
        self.blue = blue
    
    def matrices(self):
        rgb_image = np.dstack((self.red, self.green, self.blue))
        print("Displaying the pixel matrices over the Image")
        plt.imshow(rgb_image)
        plt.show()
        return self.red, self.green, self.blue
    
    def novelOperator(self,a,b,p):
        return ((a+1)*(b+1)%p)-1
    
    def MixRows(self, image, IV):
        mixedImage = image.copy()
        for i in range(0,self.m):
            mixedImage[i][1] = self.novelOperator(image[i][1],IV[i],257)
            for j in range(1,self.n):
                mixedImage[i][j] = self.novelOperator(image[i][j],mixedImage[i][j-1],257)
        return mixedImage
    
    def is_prime(self, n):
        if n <= 1:
            return False
        elif n <= 3:
            return True
        elif n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    def largest_prime_before(self, limit):
        for num in range(limit - 1, 1, -1):
            if self.is_prime(num):
                return num
            
    def ACM(self,x,y,n):
        xydash = (np.dot([[2,1],[1,1]],np.transpose([[x,y]])) % n)
        return xydash[0],xydash[1]
    
    def applyACM(self,mixedImageR,mixedImageG,mixedImageB):
        m=self.m
        n=self.n
        N = min(m,n)
        alpha = math.ceil(max(m,n)/N)
        L = N-(max(m,n) % N)
        if(alpha>1):
            eta = math.floor(L/(alpha-1))
        else:
            eta = math.floor(L/(alpha))
        k = self.largest_prime_before(int(N/2))
        x0,y0 = 0,0
        print('m= {}, n= {},\nN(length of side) = {},\nalpha(number of squares) = {},\n'
        'L(an extra length) = {},\neta(length of overlapping except last) = {}, k = {}.'
              .format(m,n,N,alpha,max(m,n)-L,eta,k))

        imgacmR = mixedImageR.copy()
        imgacmG = mixedImageG.copy()
        imgacmB = mixedImageB.copy()
        for i in range(0,alpha):
            for x in range(k):
                print(x,end=" ")
                for a in range(x0,x0+N-1):
                    for b in range(y0,y0+N-1):
                        o,p = self.ACM(a,b,N)
                        imgacmR[o[0]][p[0]],imgacmR[a][b] = imgacmR[a][b],imgacmR[o[0]][p[0]] #shuffling the values of the matrix
                        imgacmG[o[0]][p[0]],imgacmG[a][b] = imgacmG[a][b],imgacmG[o[0]][p[0]]
                        imgacmB[o[0]][p[0]],imgacmB[a][b] = imgacmB[a][b],imgacmB[o[0]][p[0]]

                if i==(alpha-1):
                    if m>n:
                        x0 = m-N
                    else:
                        y0 = n-N
                else:
                    if m > n:
                        x0 = x0 + N - eta
                    else:
                        y0 = x0 + N - eta
            return imgacmR.astype(np.uint8),imgacmG.astype(np.uint8),imgacmB.astype(np.uint8)
        
    def Mixing(self, image):
        m=self.m
        n=self.n
        for x in range(0,m):
            for y in range(1,n):
                image[x][y] = self.novelOperator(image[x][y],image[x][y-1],257)
        for y in range(0,n):
            for x in range(1,m):
                image[x][y] = self.novelOperator(image[x][y],image[x-1][y],257)
        for x in range(0,m):
            for y in range(n-2,0,-1):
                image[x][y] = self.novelOperator(image[x][y],image[x][y+1],257)
        for y in range(0,n):
            for x in range(m-2,0,-1):
                image[x][y] = self.novelOperator(image[x][y],image[x+1][y],257)

        return image.astype(np.uint8)
    
    def binary(self,num):
        bits = 8
        br = bin(num)[2:]
        # Pad the binary string with leading zeros to make it 8 bits long
        br = br.zfill(bits)

        return br
    
    def Encoding(self,image,RM):
        m=self.m
        n=self.n
        encodedImage = []
        for x in range(0,len(image)):
            count=0
            for y in range(0,len(image[0])):
                bnry = self.binary(image[x][y])
                lst=[]
                for z in range(0,len(bnry),2):
                    str=bnry[z]+bnry[z+1]
                    list = self.EncodingRules.get(str)
                    lst.append(list[(RM[(x*n+y)][count%4])-1])
                    count = count + 1
                encodedImage.append(lst)
        return encodedImage
    
    srule1=pd.DataFrame({'A':['A','C','G','T'],'C':['C','T','A','G'],'G':['G','A','T','C'],'T':['T','G','C','A']},index=['A','C','G','T'])
    srule2=pd.DataFrame({'A':['A','C','G','T'],'C':['C','T','A','G'],'G':['G','A','T','C'],'T':['T','G','C','A']},index=['A','C','G','T'])
    srule3=pd.DataFrame({'A':['G','A','T','C'],'C':['A','C','G','T'],'G':['T','G','C','A'],'T':['C','T','A','G']},index=['A','C','G','T'])
    srule4=pd.DataFrame({'A':['C','T','A','G'],'C':['T','G','C','A'],'G':['A','C','G','T'],'T':['G','A','T','C']},index=['A','C','G','T'])
    srule5=pd.DataFrame({'A':['G','A','T','C'],'C':['A','C','G','T'],'G':['T','G','C','A'],'T':['C','T','A','G']},index=['A','C','G','T'])
    srule6=pd.DataFrame({'A':['C','T','A','G'],'C':['T','G','C','A'],'G':['A','C','G','T'],'T':['G','A','T','C']},index=['A','C','G','T'])
    srule7=pd.DataFrame({'A':['T','G','C','A'],'C':['G','A','T','C'],'G':['C','T','A','G'],'T':['A','C','G','T']},index=['A','C','G','T'])
    srule8=pd.DataFrame({'A':['T','G','C','A'],'C':['G','A','T','C'],'G':['C','T','A','G'],'T':['A','C','G','T']},index=['A','C','G','T'])

    def getRule(self,rule,dM,image):
        if(rule==0):
            return self.srule1[dM][image]
        if(rule==1):
            return self.srule2[dM][image]
        if(rule==2):
            return self.srule3[dM][image]
        if(rule==3):
            return self.srule4[dM][image]
        if(rule==4):
            return self.srule5[dM][image]
        if(rule==5):
            return self.srule6[dM][image]
        if(rule==6):
            return self.srule7[dM][image]
        if(rule==7):
            return self.srule8[dM][image]
        
    def Substitution(self,encodedImage,DSM,RM):
        substituted = []
        for i in range(len(encodedImage)):
            lst=[]
            for j in range(len(encodedImage[0])):
                dM = DSM[(i*len(encodedImage[0]))+j]
                rule = RM[i][j]-1
                dI = self.getRule(rule,dM,encodedImage[i][j])
                lst.append(dI)
            substituted.append(lst)
        return substituted
    
    rule1={'A':'00','C':'01','G':'10','T':'11'}
    rule2={'A':'00','G':'01','C':'10','T':'11'}
    rule3={'C':'00','A':'01','T':'10','G':'11'}
    rule4={'G':'00','A':'01','T':'10','C':'11'}
    rule5={'C':'00','T':'01','A':'10','G':'11'}
    rule6={'G':'00','T':'01','A':'10','C':'11'}
    rule7={'T':'00','C':'01','G':'10','A':'11'}
    rule8={'T':'00','G':'01','C':'10','A':'11'}
    
    def getBinEqui(self,rul,dD):
        if(rul==0):
            return self.rule1[dD]
        if(rul==1):
            return self.rule2[dD]
        if(rul==2):
            return self.rule3[dD]
        if(rul==3):
            return self.rule4[dD]
        if(rul==4):
            return self.rule5[dD]
        if(rul==5):
            return self.rule6[dD]
        if(rul==6):
            return self.rule7[dD]
        if(rul==7):
            return self.rule8[dD]

    def getDecEqui(self,n):
        return int(n,2)
    
    def Decoding(self,image,RM):
        Pfinal = []
        for i in range(len(image)):
            de = 0
            for j in range(4):
                dD = image[i][j]
                rules = RM[i][j]-1
                u = self.getBinEqui(rules,dD)
                u = self.getDecEqui(u)
                de = 4*de+u
            Pfinal.append(de)
        return Pfinal

    def toMatrix(self,lst):
        m=self.m
        n=self.n
        final = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                final[i][j] = lst[(i*n)+j]
        return final.astype(np.uint8)
    
    def getFinalCipherImage(self,finalR,finalG,finalB):
        rgb_image = np.dstack((finalB,finalG,finalR,))
        plt.imshow(rgb_image)
        plt.show()
        return cv2.merge((finalB,finalG,finalR))
    
    def mixRows(self,r,g,b,IV):
        mixRed = self.MixRows(r,IV)
        mixGreen = self.MixRows(g,IV)
        mixBlue = self.MixRows(b,IV)
        return mixRed,mixGreen,mixBlue
    
    def mixing(self,r,g,b):
        mixR = self.Mixing(r)
        mixG = self.Mixing(g)
        mixB = self.Mixing(b)
        return mixR,mixG,mixB
    
    def encoding(self,r,g,b,RM):
        ecdR = self.Encoding(r,RM)
        ecdG = self.Encoding(g,RM)
        ecdB = self.Encoding(b,RM)
        return ecdR,ecdG,ecdB
    
    def substitution(self,r,g,b,DSM,RM):
        subR = self.Substitution(r,DSM,RM)
        subG = self.Substitution(g,DSM,RM)
        subB = self.Substitution(b,DSM,RM)
        return subR,subG,subB
    
    def decoding(self,r,g,b,RM):
        decR = self.Decoding(r,RM)
        decG = self.Decoding(g,RM)
        decB = self.Decoding(b,RM)
        return decR,decG,decB
    
    def applyEncryption(self,image,IV,RM,DSM):
        b,g,r = cv2.split(image)
        print("Applying Algorithm Steps one by one\n")
        print("Steps 1: mix rows\n")
        mrR,mrG,mrB = self.mixRows(r,g,b,IV)
        print("Step 2: applying ACM algorithm\n")
        acmR,acmG,acmB = self.applyACM(mrR,mrG,mrB)
        print("\nStep 3: applying mixing \n")
        mixR,mixG,mixB = self.mixing(acmR,acmG,acmB)
        print("Step 4: applying Encoding\n")
        encR,encG,encB = self.encoding(mixR,mixG,mixB,RM)
        print("Step 5: applying Substitution\n")
        subR,subG,subB = self.substitution(encR,encG,encB,DSM,RM)
        print("Step 6: applying Decoding")
        decR,decG,decB = self.decoding(subR,subG,subB,RM)
        finalR = self.toMatrix(decR)
        finalG = self.toMatrix(decG)
        finalB = self.toMatrix(decB)
        final_cipher_image = self.getFinalCipherImage(finalR,finalG,finalB)
        return final_cipher_image
        
        
        


# In[13]:


class Decryption(Encryption):
    def __init__(self, cipherImage):
        self.m = cipherImage.shape[0]
        self.n = cipherImage.shape[1]
        blue, green, red = cv2.split(cipherImage)
        self.red = red
        self.green = green
        self.blue = blue
        
        
    revrule1=pd.DataFrame({'A':['A','C','G','T'],'C':['G','A','T','C'],'G':['C','T','A','G'],'T':['T','G','C','A']},index=['A','C','G','T'])
    revrule2=pd.DataFrame({'A':['A','C','G','T'],'C':['G','A','T','C'],'G':['C','T','A','G'],'T':['T','G','C','A']},index=['A','C','G','T'])
    revrule3=pd.DataFrame({'A':['C','T','A','G'],'C':['A','C','G','T'],'G':['T','G','C','A'],'T':['G','A','T','C']},index=['A','C','G','T'])
    revrule4=pd.DataFrame({'A':['G','A','T','C'],'C':['T','G','C','A'],'G':['A','C','G','T'],'T':['C','T','A','G']},index=['A','C','G','T'])
    revrule5=pd.DataFrame({'A':['C','T','A','G'],'C':['A','C','G','T'],'G':['T','G','C','A'],'T':['G','A','T','C']},index=['A','C','G','T'])
    revrule6=pd.DataFrame({'A':['G','A','T','C'],'C':['T','G','C','A'],'G':['A','C','G','T'],'T':['C','T','A','G']},index=['A','C','G','T'])
    revrule7=pd.DataFrame({'A':['T','G','C','A'],'C':['C','T','A','G'],'G':['G','A','T','C'],'T':['A','C','G','T']},index=['A','C','G','T'])
    revrule8=pd.DataFrame({'A':['T','G','C','A'],'C':['C','T','A','G'],'G':['G','A','T','C'],'T':['A','C','G','T']},index=['A','C','G','T'])
    
    def getRevRule(self, rule,dM,image):
        if(rule==0):
            return self.revrule1[dM][image]
        if(rule==1):
            return self.revrule2[dM][image]
        if(rule==2):
            return self.revrule3[dM][image]
        if(rule==3):
            return self.revrule4[dM][image]
        if(rule==4):
            return self.revrule5[dM][image]
        if(rule==5):
            return self.revrule6[dM][image]
        if(rule==6):
            return self.revrule7[dM][image]
        if(rule==7):
            return self.revrule8[dM][image]
        
    def rSubstitution(self,encodedImage,DSM,RM):
        substituted = []
        for i in range(len(encodedImage)):
            lst=[]
            for j in range(len(encodedImage[0])):
                dM = DSM[(i*len(encodedImage[0]))+j]
                rule = RM[i][j]
                dI = self.getRevRule(rule-1,dM,encodedImage[i][j])
                lst.append(dI)
            substituted.append(lst)
        return substituted
    
    def reverseNovelOperator(self,result, bb, par):
        ans = ((result + 1) * self.modinv(bb + 1, par)) % par - 1
        return ans

    # Modular multiplicative inverse function
    def modinv(self,aa, t):
        m0, x0, x1 = t, 0, 1
        while aa > 1:
            q = aa // t
            t, aa = aa % t, t
            x0, x1 = x1 - q * x0, x0
        return x1 + m0 if x1 < 0 else x1

    def rMixing(self,image):
        m=self.m
        n=self.n
        for y in range(n-1,-1,-1):
            for x in range(1,m-2):
                image[x][y] = self.reverseNovelOperator(image[x][y],image[x+1][y],257)
        for x in range(m-1,-1,-1):
            for y in range(1,n-2):
                image[x][y] = self.reverseNovelOperator(image[x][y],image[x][y+1],257)
        for y in range(n-1,-1,-1):
            for x in range(m-1,0,-1):
                image[x][y] = self.reverseNovelOperator(image[x][y],image[x-1][y],257)
        for x in range(m-1,-1,-1):
            for y in range(n-1,0,-1):
                image[x][y] = self.reverseNovelOperator(image[x][y],image[x][y-1],257)

        return image.astype(np.uint8)


    def applyRevACM(self,ImageR,ImageG,ImageB):
        m=self.m
        n=self.n
        N = min(m,n)
        alpha = math.ceil(max(m,n)/N)
        L = N-(max(m,n) % N)
        if(alpha>1):
            eta = math.floor(L/(alpha-1))
        else:
            eta = math.floor(L/(alpha))
        k = self.largest_prime_before(int(N/2))
        x0,y0 = 0,0
        print('m= {}, n= {},\nN(length of side) = {},\nalpha(number of squares) = {},\n'
        'L(an extra length) = {},\neta(length of overlapping except last) = {}, k = {}.'
              .format(m,n,N,alpha,max(m,n)-L,eta,k))

        imgacmR = ImageR.copy()
        imgacmG = ImageG.copy()
        imgacmB = ImageB.copy()
        for i in range(0,alpha):
            for x in range(k):
                print(x,end=" ")
                for a in range(x0 + N-2, x0 - 1, -1):
                    for b in range(y0 + N-2, y0 - 1, -1):
                        o,p = self.ACM(a,b,N)
                        imgacmR[o[0]][p[0]],imgacmR[a][b] = imgacmR[a][b],imgacmR[o[0]][p[0]] #shuffling the values of the matrix
                        imgacmG[o[0]][p[0]],imgacmG[a][b] = imgacmG[a][b],imgacmG[o[0]][p[0]]
                        imgacmB[o[0]][p[0]],imgacmB[a][b] = imgacmB[a][b],imgacmB[o[0]][p[0]]

                if i==(alpha-1):
                    if m>n:
                        x0 = m-N
                    else:
                        y0 = n-N
                else:
                    if m > n:
                        x0 = x0 + N - eta
                    else:
                        y0 = x0 + N - eta
            return imgacmR.astype(np.uint8),imgacmG.astype(np.uint8),imgacmB.astype(np.uint8)


    def rmixRow(self,image,IV):
        # Iterate over the rows in reverse order
        m=self.m
        n=self.n
        Dimg = image.copy()
        for i in range(m - 1, -1, -1):
            for j in range(n - 1, 1, -1):
                Dimg[i][j] = self.reverseNovelOperator(image[i][j], image[i][j-1], 257)
            Dimg[i][1] = self.reverseNovelOperator(image[i][1], IV[i], 257)

        # The second column
        for j in range(n - 1, -1, -1):
            Dimg[j][1] = Dimg[j][0]
        return Dimg.astype(np.uint8)
    
    def getDecodedImage(self,dred,dgreen,dblue):
        super().getFinalCipherImage(dred,dgreen,dblue)
        
    
    def revSubstitution(self,r,g,b,DSM,RM):
        rsR = self.rSubstitution(r,DSM,RM)
        rsG = self.rSubstitution(g,DSM,RM)
        rsB = self.rSubstitution(b,DSM,RM)
        return rsR,rsG,rsB
    
    def revMixing(self,r,g,b):
        rmR = self.rMixing(r)
        rmG = self.rMixing(g)
        rmB = self.rMixing(b)
        return rmR,rmG,rmB
    
    def revmixRows(self,r,g,b,IV):
        rmrR = self.rmixRow(r,IV)
        rmrG = self.rmixRow(g,IV)
        rmrB = self.rmixRow(b,IV)
        return rmrR,rmrG,rmrB
    
    def applyDecryption(self,image,IV,RM,DSM):
        print("Applying Decryption Algorithm on provided image\n")
        b,g,r = cv2.split(image)
        print("Step 1: applying Encoding\n")
        enR,enG,enB = self.encoding(r,g,b,RM)
        print("Step 2: applying Reverse Substitution\n")
        rsR,rsG,rsB = self.revSubstitution(enR,enG,enB,DSM,RM)
        print("Step 3: applying Decoding\n")
        dcR,dcG,dcB = self.decoding(rsR,rsG,rsB,RM)
        fred = self.toMatrix(dcR)
        fgreen = self.toMatrix(dcG)
        fblue = self.toMatrix(dcB)
        print("Step 4: applying reverse mixing\n")
        mixR,mixG,mixB = self.revMixing(fred,fgreen,fblue)
        print("Step 5: applying reverse ACM algorithm\n")
        acmR,acmG,acmB = self.applyRevACM(mixR,mixG,mixB)
        print("\nStep 6: applying reverse mix rows\n")
        mrR,mrG,mrB = self.revmixRows(acmR,acmG,acmB,IV)
        decodedImage = self.getDecodedImage(mrR,mrG,mrB)
        return decodedImage
        
    


# In[ ]:




