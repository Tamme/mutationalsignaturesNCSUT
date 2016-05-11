#! /usr/bin/python

# usage: python topics.py <beta file> <vocab file> <num words>
#
# <beta file> is output from the lda-c code
# <vocab file> is a list of words, one per line
# <num words> is the number of words to print from each topic

import sys
from math import exp


def print_topics(beta_file, vocab_file, nwords = 25):

    # get the vocabulary

    vocab = file(vocab_file, 'r').readlines()
    # vocab = map(lambda x: x.split()[0], vocab)
    vocab = map(lambda x: x.strip(), vocab)
    #print vocab
    # for each line in the beta file

    indices = range(len(vocab))
    topic_no = 0
    fullString = ''
    for topic in file(beta_file, 'r'):
#        print 'topic %03d' % topic_no
        topic = map(float, topic.split())
#	print(len(topic))
#	return
	topicPow = [ exp(x) for x in topic ]
	strTopicPow = list(map(str, topicPow))
	fullString += ' '.join(strTopicPow) + '\n'
#	print string
#	return
#	print sum(topicPow2)
#	return
	#print(indices)
        indices.sort(lambda x,y: -cmp(topic[x], topic[y]))
        for i in range(nwords):
            print '   %s' % vocab[indices[i]], '   %s' % exp(topic[indices[i]])
        topic_no = topic_no + 1
        print '\n'
    fil = open(beta_file + '.prob.txt', 'w')
    fil.write(fullString)
    fil.close()
    print "made new"

if (__name__ == '__main__'):

    if (len(sys.argv) != 4):
       print 'usage: python topics.py <beta-file> <vocab-file> <num words>\n'
       sys.exit(1)

    beta_file = sys.argv[1]
    vocab_file = sys.argv[2]
    nwords = int(sys.argv[3])
    print_topics(beta_file, vocab_file, nwords)
