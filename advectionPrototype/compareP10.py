from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

def tiffToArray(tiffPath):
    im = Image.open(tiffPath)
    return np.array(im)

if __name__=="__main__":

    '''
    file = "I:/work-he/apps/safi/data/IBI/P10with5pctDeficit3m"
    yearDeficit = []
    with open(file + ".txt") as f:
        fich = f.readline()
        for str in fich.split('['):
            if len(str.split(']')) > 1:
                ligne = str.split(']')[0]
                yearDeficit += [[float(val) for val in ligne.split(',')]]
    yearDeficit = np.array(yearDeficit)
    '''
    withoutDeficitTiffPath = "I:/work-he/apps/safi/data/IBI/P102020.tiff"
    withDeficitTiffPath10 = "I:/work-he/apps/safi/data/IBI/P10with10pctDeficit3m.tiff"
    withDeficitTiffPath5 = "I:/work-he/apps/safi/data/IBI/P10with5pctDeficit3m.tiff"
    withDeficitTiffPath2 = "I:/work-he/apps/safi/data/IBI/P10with2pctDeficit3m.tiff"

    NoDeficitP10 = tiffToArray(withoutDeficitTiffPath)
    #newDecile = np.sum((yearDeficit>NoDeficitP10)*1.,axis=0)/len(yearDeficit)

    DeficitArray10 = tiffToArray(withDeficitTiffPath10)
    DeficitArray5 = tiffToArray(withDeficitTiffPath5)
    DeficitArray2 = tiffToArray(withDeficitTiffPath2)

    #diff5 = DeficitArray5-DeficitArray10
    diff10 = NoDeficitP10-DeficitArray10
    diff5 = NoDeficitP10 - DeficitArray5
    diff2 = NoDeficitP10 - DeficitArray2

    fig4, ax4 = plt.subplots()
    plt.imshow(NoDeficitP10)
    ax4.set_title('NO3 P10 without deficit')
    plt.clim(1, 8)
    plt.colorbar()
    plt.show()

    fig5, ax5 = plt.subplots()
    plt.imshow(DeficitArray10)
    ax5.set_title('NO3 P10 with 10% deficit')
    plt.clim(1, 8)
    plt.colorbar()
    plt.show()

    fig6, ax6 = plt.subplots()
    plt.imshow(DeficitArray5)
    ax6.set_title('NO3 P10 with 5% deficit')
    plt.clim(1, 8)
    plt.colorbar()
    plt.show()

    fig7, ax7 = plt.subplots()
    plt.imshow(DeficitArray2)
    ax7.set_title('NO3 P10 with 2% deficit')
    plt.clim(1, 8)
    plt.colorbar()
    plt.show()

    fig1, ax1 = plt.subplots()
    plt.imshow(diff10)
    ax1.set_title('Difference between NO3 P10 with 10% deficit and NO3 P10 without deficit')
    plt.clim(0, 1.5)
    plt.colorbar()
    plt.show()

    fig2, ax2 = plt.subplots()
    plt.imshow(diff5)
    ax2.set_title('Difference between NO3 P10 with 5% deficit and NO3 P10 without deficit')
    plt.clim(0, 1.5)
    plt.colorbar()
    plt.show()

    fig3, ax3 = plt.subplots()
    plt.imshow(diff2)
    ax3.set_title('Difference between NO3 P10 with 2% deficit and NO3 P10 without deficit')
    plt.clim(0, 1.5)
    plt.colorbar()
    plt.show()

    maskposition = np.where(np.isnan(DeficitArray10))
    negval = (DeficitArray10 < 0) * -1.
    negval[maskposition] = np.nan
    negval = negval * DeficitArray10
    fig, ax = plt.subplots()
    plt.imshow(negval)
    ax.set_title('NO3 P10 10% deficit IBI, negatives values')
    plt.colorbar()
    plt.clim(0, 0.05)
    plt.show()

    maskposition = np.where(np.isnan(DeficitArray10))
    negval = (DeficitArray5 < 0) * -1.
    negval[maskposition] = np.nan
    negval = negval * DeficitArray5
    fig, ax = plt.subplots()
    plt.imshow(negval)
    ax.set_title('NO3 P10 5% deficit IBI, negatives values')
    plt.colorbar()
    plt.clim(0, 0.05)
    plt.show()

    maskposition = np.where(np.isnan(DeficitArray10))
    negval = (DeficitArray2 < 0) * -1.
    negval[maskposition] = np.nan
    negval = negval * DeficitArray2
    fig, ax = plt.subplots()
    plt.imshow(negval)
    ax.set_title('NO3 P10 2% deficit IBI, negatives values')
    plt.colorbar()
    plt.clim(0, 0.05)
    plt.show()

    '''
    fig1, ax2 = plt.subplots()
    plt.imshow(diff10)
    ax2.set_title('Difference between with 10% and without deficit')
    plt.clim(0, 0.2)
    plt.colorbar()
    plt.show()'''