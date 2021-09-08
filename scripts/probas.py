import numpy
import matplotlib.pyplot as plt

def pb_overlap(n, A, B):
    if A + B > n: return 1
    return 1 - numpy.prod([(n-A-i)/(n-i) for i in range(B)])

plt.figure(figsize=(12, 3))
As = [1, 2, 3, 4, 5, 7, 10, 15, 20, 25]
ns = [25, 50, 100]
for ni, n in enumerate(ns):
    plt.subplot(1, len(ns), ni+1)
    B = numpy.arange(1, 26)
    for A in As:
        ps = numpy.array([pb_overlap(n, A, s) for s in B])
        plt.plot(B, ps, '.-', label="A = %d" % A)    

    plt.xticks([1, 5, 10, 15, 20, 25])
    plt.ylim([-.01, 1.01])
    plt.xlabel("B")
    if ni == 0:
        plt.ylabel("Overlap probability")
    plt.title("n = %d" % n)
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.subplots_adjust(right=0.8, bottom=0.2, left=0.1, top=0.9)
plt.show()
