# project2

Βαμβουρέλλης Ευστράτιος 1115201600014
Αποστόλου Φίλιππος 1115201600007

--------------------------------ΠΡΟΣΟΧΗ---------------------------------------
1: Για την ανάγνωση αρχείων χρησιμοποιήσαμε την παλιά συνάρτηση από την πρώτη εργασία. Τα αρχεία input πρέπει να είναι σε λίγο διαφορετικό φορμάτ από αυτό που ζητάει η εκφώνηση!! Φτιάξαμε ένα σκριπτ που μετατρέπει ένα αρχείο στο ζητούμενο φορμάτ.
2: Το όρισμα -complete στο δικό μας πρόγραμμα ονομάζεται -a 0 η 1
3: Στο αρχείο .conf πρέπει να δοθεί flag που ορίζει αν το input είναι καμπύλες η διανύσματα.
4: Στο .out αρχείο εμφανίζονται κάποια κόμματα στο τέλος της εκτύπωσης των vector η των curve η του sihlouette. Αυτά ΔΕΝ είναι λάθος εκτέλεση του προγράμματος (απλά δεν δουλεύει το \b).
5: Ένα cluster μπορεί να μην έχει κανένα σημείο. Αυτό συμβαίνει μόνο όταν γίνεται χρίση του update2 γιατί όταν κατα τη διάρκεια του clastering ένα cluster δεν έχει κανένα σημείο ΔΕΝ αλλάζουμε το κέντρο του (πχ με ένα άλλο τυχαίο σημείο). Αφήνουμε το κέντρο το ίδιο και στη συνέχεια μπορεί το cluster αυτό να αποκτήσει και πάλι σημεία. Υπάρχει όμως πάντα η πιθανότητα να τερματίσει άδειο.
6: Στην παραπάνω περίπτωση το sihlouette επιστρέφει την τιμή -2 ενδεικτικά για να μας ενημερώσει ότι το cluster είναι κενό (το -2 πιάνεται στον μέσο όρο οπότε μπορεί να τον αλλοιώσει αρκετά).

--------------------------------.conf---------------------------------------

number_of_clusters: unsigned int αριθμός cluster k
number_of_grids: unsigned int αριθμός διαφορετικών grid (l for lsh_curves)
number_of_vector_hash_tables: unsigned int l for lsh_vectors
number_of_vector_hash_functions: unsigned int k for lsh (nomber of f’s in g)
lsh_window: double >0 window for lsh
lsh_curves_grids_delta: double >0 η εμτατόπιση δ στα grids του lsh_curves
input_contains_vectors: 0 η 1 αναλόγως αν το input περιέχει vectors η curves
lsh_multimap_container_size: unsigend int lsh container size
stop_when_centers_dont_change: 0 η 1 ανάλογα με την συνθήκη τερματισμού που επιλέγουμε
curve_tolerance: 0.001
vector_tolerance: 0.01
center_tolerance:
max_iterations: unsigned int ο μέγιστος αριθμός επαναλήψεων που θα εκτελέσει ο αλγόριθμος
lsh_curves_pad_number: double lsh_curves pad number
max_curve_size: τα μέγιστο μήκος καμπύλης των δεδομένων
brute_update_1: 0 η 1 αναλόγως με την παραλλαγή του update1 που επιλέγουμε

---------------------------------παραδωταία--------------------------------------

1: Υλοποιήσαμε όλους τους αλγόριθμους της εκφώνησης και κάποιους έξτρα με δικές μας παραλλαγές και βελτιώσεις. Όλοι εκτελούνται από μια main με την σειρά.

Μεταγλώττιση: make
Εκτέλεση: ./main – i input file -c file.conf -o file.out -a 0

Το -a 1 αντιστοιχεί στο -complete της εκφώνησης. Το -a 0 αντιστοιχεί στην απουσία -complete της εκφώνησης (μπορεί να παραληφθεί καθώς έχει default τιμή).

2: unit testing για 3 συναρτήσεις βλ παρακάτω.

3: Ένα σκριπτ format input. Αυτό το σκριπτ δουλεύει μόνο για τα αρχεία δεδομένων που είναι ανεβασμένα στο eclass. Μετατρέπει το αρχείο στο φορματ των αρχείων της πρώτης εργασίας, για να μπορούμε να τα διαβάσουμε χρησιμοποιώντας την ίδια συνάρτηση. ΠΡΟΣΟΧΗ δεν έχει δοκιμαστεί για άλλα αρχεία input (πχ αν ο χωρισμός δεν γίνεται με \t αλλά με ,)

4: makefile το οποίο κάνει χρήση των flags -ggdb -g -m64 -O0.

--------------------------------update1_brute---------------------------------------

Έχουμε υλοποιήσει 2 παραλλαγές του update1. Υπάρχει η επιλογή update στο .conf.
1: update1_brute είναι ο αλγόριθμος που προτείνεται απο την εκφώνηση στην brutal force υλοποίηση του.
2: update1 είναι μια βελτίωση στον αλγόριθμο brute χρησιμοποιώντας lsh. Συνοπτικά βρίσκουμε τον μέσο όρο των σημείων (όπως στο update2)  και μετά επισημοποιώντας lsh βρίσκουμε το πιο κοντινό “πραγματικό” σημείο (σημείο απο το dataset). Είναι μαθηματικά προφανές οτι ο μέσος όρος των σημείων είναι αυτός που ελαχιστοποιεί το άθροισμα των αποστάσεων, οπότε βρίσκοντας το πιο κοντινό σημείο στον μέσο όρο καταλήγουμε στο ίδιο αποτέλεσμα που θα είχαμε με το brute force, αλλά σε πολύ καλύτερο χρόνο. (δεν έχουμε πάντα το ίδιο αποτέλεσμα, στην περίπτωση των καμπυλών βλ σχόλια)

--------------------------------συνθήκες τερματισμού---------------------------------------

Έχουμε υλοποιήσει δύο συνθήκες τερματισμού. Δίνεται η δυνατότητα επιλογής  μέσο του stop_when_centers_dont_change. Επιπλέον, υπάρχουν και τα max_iterations όπου ο αλγόριθμος σταματάει είτε έχει φτάσει σε σύγκλιση είτε οχι.
1: centers_dont_change ελέγχει αν τα νέα κέντρα είναι “κοντά” στα παλιά. Το “κοντά” ορίζεται απο το center_ tolerance.
2: clusters_dont_change ελέγχει αν τα νέα cluster είναι ίδια με τα παλιά.


--------------------------------initialization1---------------------------------------

Έχουμε υλοποιήσει 2 τρόπους για την επιλογή των τυχαίων κέντρων. Ο ένας είναι ο “μαθηματικά σωστότερος” και ο άλλος ο “προγραμματιστικά βέλτιστος”. Ο πρώτος είναι σε σχόλια. Ο πρώτος δεν είναι αποδεκτός σε μια πραγματική εφαρμογή γιατί σπαταλάει πολύ περισσότερη μνήμη και χρόνο απο τον δεύτερο. Ο δεύτερος σε πραγματικές συνθήκες όπου το k είναι πολύ μικρότερο από το πλήθος των δεδομένων είναι πολύ καλύτερος.


--------------------------------curve_get_mean---------------------------------------

 Στην αρχικοποίηση του C έχουμε υλοποιήσει δυο τρόπους. Αυτόν που προτείνει η εκφώνηση και έναν δικό μας λίγο διαφορετικό. Ο δικός μας είναι σε σχόλια. Ο δικός μας τρόπος ουσιαστικά δεν ψάχνει μια τυχαία καμπύλη με μήκος >=λ, αλλά παίρνει την πρώτη καμπύλη που βρίσκει, με μήκος >=λ, ψάχνοντας γραμμικά. Παρατηρήσαμε, μετά απο δωκιμές, οτι με αυτή την παραλλαγή έχουμε λίγο πιο γρήγορη εύρεση του μέσου όρου και χρησιμοποιώντας λίγο λιγότερη μνήμη. Η επίπτωση που έχει αυτό στην δημιουργία των τελικών clster είναι περίπου 0.1-0.25 με την μετρική sihloutte.

---------------------------------initialization2--------------------------------------

Στο initialization2 υπήρχε μια μικρή ασάφεια στις διαφάνειες, σχετικά με τον αριθμό χ και το διάστημα στο ποίο πρέπει να τον επιλέξουμε. Έχουμε υλοποιήσει και τους δύο δυνατούς τρόπους και σε σχόλια έχουμε τον 1.
1: πρέπει το διάστημα να περιέχει το 0 (η αριστερή ανισότητα των διαφανειών είναι λάθος).
2: το διάστημα δεν πρέπει να περιέχει το 0. (το διάστημα επιλογής του χ στις διαφάνειες είναι λάθος)
--------------------------------Unit_testing ---------------------------------------
Χρησιμοποιήσαμε το unit testing της google παραθέτω τις εντολές εγκατάστασης:

sudo apt-get install libgtest-dev
sudo apt-get install cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
sudo cp *.a /usr/lib
Για μεταγλώτηση στο directory που  έχουμε το unit_testing_main.cpp  

cmake CMakeLists.txt
make
./executeTests
Τα παραδειγματα που τρεχουμε ειναι τα εξης :
1.Το binary search για εναν αξοντα πινακα διαστηματων για τις περιπτωσεις που του στοιχειο ειναι αριστερα ,δεξια του κεντρου καθως κι τις οριακες περιπτωσεις τερμα δεξια κι αριστερα.
2.Το mean vectors βρισκει το μεσο vector κι τα παραδειγματα ειναι για δυο vector μηκους 10 κι για πολλαπλά vector μηκους 10 καποια απο αυτα να ειναι ισα με το μηδεν κι ολα τα αλλα ειναι το δυανισμα [0,1,2,3,4,5,6,7,8,9](μπορει να αλλαχθει ο αριθμος τους μεσα πο την main)
3.Το assigment2 το οποίο  το τσεκάρουμε για κ=2 σε ενα σετ  που εχει πολλα σημεια στο κεντρο των αξονων κι 2 απομακρυσμένα σημεία . Η μεγίστη απόσταση μανχαταν μεταξυ των σημείων κοντα στο κεντρο ειναι 4 . Επειδη ο αλγοριθμος ειναι τυχαιος ο ελενχος μας δεν είναι αυστηρός δλδ καποιες φορες δεν βρισκει κάποιο από αυτά σαν κέντρο.

---------------------------------σχόλια--------------------------------------

1: Όλοι οι αλγόριθμοι έχουν πολύ καλά αποτελέσματα στα vectors. Με ένα μέσο όρο 0.97.
2: Το initialasation2 είναι ξεκάθαρα πολύ ανώτερο από το initialasation1.
3: Στα curves έχουμε λίγο χειρότερα αποτελέσματα, περίπου μ.ο. 0.4
4: Οι αλγόριθμοι για vectors είναι λίγο ευαίσθητοι τις παραμέτρους σύγκλισης (tolerance). Τα curves είναι πολύ ευαίσθητα σε αυτές τις παραμέτρους.
5: Είναι δύσκολο να επιλέξεις καλές παραμέτρους σύγκλισης για τα curves.
6: Το stop_when_equal_clusters έχει πολύ καλύτερα αποτελέσματα, αν και είναι πολύ πιο αργό, για τα curves.
7: Η δικιά μας παραλλαγή του update1 έχει σημαντικά χειρότερα αποτελέσματα στα curves (μέχρι και 0.4 διαφορά). Αυτό είναι λογικό καθώς το μέσο curve δεν είναι “ακριβώς μαθηματικά ορισμένο”.
8: Τα curves έχουν πολύ χειρότερα αποτελέσματα με το update2. Αυτό είναι λογικό για τον λόγο που προανέφερετε.
9: Με τις σωστές παραμέτρους το assigment2 είναι πολύ πιο γρήγορο απο το assigment1.
