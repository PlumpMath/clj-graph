(ns clj-graph.core
  (:require [clojure.java.io :refer [reader]]
            [clojure.set :as set]
            [clojure.core.matrix.stats :as s]))

(defn read-graph [file]
  (with-open [rdr (reader file)]
    (->>
     (for [line (line-seq rdr)
           :let [[[_ a b]] (re-seq #"(\d+) (\d+)" line)]
           :when (and a b)
           :let [a (Integer/parseInt a)
                 b (Integer/parseInt b)]]
       [a b])
     (reduce (fn [old [a b]]
               (-> old
                   (update-in [a] #(conj (or % (sorted-set)) b))
                   (update-in [b] #(conj (or % (sorted-set)) a))))
             {}))))

(def g (read-graph "Florida-bay.txt"))

;; Test graph from wikipedia
(def tg {1 #{3}
         2 #{3}
         3 #{1 2 4 5}
         4 #{3 5}
         5 #{3 4}})

(defn degree-sequence [G v-sub]
  (let [s (select-keys G v-sub)]
    (->> (for [[v ns] s]
           (count (filter v-sub ns)))
         sort
         vec)))

(degree-sequence tg #{1 2 3})

;; d)

(defn enumerate-subgraphs
  "Enumeration of ESU (FANMOD)

  This enumerates all connected subgraphs of G of size k with the ESU
  algorithm:

  Wernicke S (2006). \"Efficient detection of network motifs\". IEEE/ACM
  Transactions on Computational Biology and Bioinformatics. 3 (4): 347–359.
  doi:10.1109/tcbb.2006.51

  It is a faithful direct implementation of the pseudo code given in Wernicke
  2006 as also quoted by https://en.wikipedia.org/wiki/Network_motif

  "
  [G k]
  (loop [res []
         ;; explicit extend-subgraph stack on heap
         stack (for [v (keys G)
                     :let [v-ext (set (filter #(> % v) (G v)))]]
                 [#{v} v-ext v])]
    (if (empty? stack)
      res
      (let [[[v-sub v-ext v] & stack] stack]
        (if (= (count v-sub) k)
          (recur (conj res v-sub) stack)
          ;; else still work to do => extend stack
          (let [calls (loop [extend-subgraphs []
                             v-ext v-ext]
                        (let [w (rand-nth (seq v-ext))]
                          (if-not w
                            extend-subgraphs
                            (let [v-sub' (conj v-sub w)
                                  ;; as defined in section 2.1 of Wernicke 2006
                                  ns-vsub (into v-sub (mapcat G) v-sub)
                                  v-ext' (into (disj v-ext w)
                                               (filter #(and (> % v)
                                                             (not (ns-vsub %))))
                                               (G w))]
                              (recur (conj extend-subgraphs [v-sub' v-ext' v])
                                     (disj v-ext w))))))]
            (recur res (concat calls stack))))))))

;; G-Trie would also have looked reasonable, but way beyond the scope of this
;; exercise, ESU FANMOD was the only reasonable algorithm (performance vs.
;; effort) for this problem.


(enumerate-subgraphs tg 3)

(time
 (def g-subgraphs (enumerate-subgraphs g 4)))

(defn degree-seq-counts [subgraphs]
  (->> subgraphs
       (map (partial degree-sequence g))
       frequencies))

(def g-degree-counts (degree-seq-counts g-subgraphs))

(comment
   ;; =>
  {[1 1 2 2] 666566,
   [1 1 1 3] 529287,
   [2 2 2 2] 131661,
   [1 2 2 3] 414426,
   [2 2 3 3] 144942,
   [3 3 3 3] 14126}


  (count (enumerate-subgraphs g 4)))


;; calculate e)

(defn erdos-renyi [N m]
  (->> (for [a (range N)
             b (range N)
             :when (> b a)]
         [a b])
       shuffle
       (take m)
       (reduce (fn [g [a b]]
                 (-> g
                     (update a #(conj (or % #{}) b))
                     (update b #(conj (or % #{}) a))))
               {})))


(defn sample-node-by-degree [g]
  (let [ds (map (fn [[_ ns]] (count ns)) g)
        Z (reduce + ds)
        r (rand)]
    (loop [[[n nghbs] & ns] (seq g)
           s 0.0] ;; accumulate
      (let [pn (/ (count nghbs)
                  Z)]
        (if (> (+ s pn) r)
          n
          (recur ns (+ s pn)))))))


(defn add-node [g [n d]]
  (let [nghbs (->> (repeatedly #(sample-node-by-degree g))
                   distinct
                   (take d))]
    (-> (reduce (fn [g n']
                  (update g n' #(conj (or % #{}) n)))
                g
                nghbs)
        (assoc n (set nghbs)))))


(defn barabasi-albert [core t m]
  (loop [g core
         t t
         n (inc (apply max (keys core)))]
    (if (zero? t)
      g
      (recur (add-node g [n m])
             (dec t)
             (inc n)))))


(defn rewire-config-model [g]
  (let [ds (->> (map (fn [[n nghbs]]
                       [(count nghbs) n])
                     g)
                (reduce (fn [ds [c n]]
                          (update ds c #(conj (or % #{}) n)))
                        {}))
        all-degrees (into [] (keys ds))
        [a c b d]
        (first (for [[a b] (repeatedly #(let [dnds (seq (ds (rand-nth all-degrees)))]
                                          [(rand-nth dnds) (rand-nth dnds)]))
                     :when (and a b (not= a b))
                     :let [na (g a)
                           nb (g b)
                           c (rand-nth (seq na))
                           d (rand-nth (seq nb))]
                     :when (and c d (not= c d)
                                (not (na d))
                                (not (nb c)))]
                 [a c b d]))]
    (-> g
        (update a disj c)
        (update a conj d)
        (update b disj d)
        (update b conj c))))



(defn configuration-model [g steps]
  (reduce (fn [g _] (rewire-config-model g))
          g
          (range steps)))



(comment
  (rewire-config-model {1 #{2 3}
                        2 #{1 5}
                        3 #{1 2 4}
                        4 #{3}
                        5 #{2}})

  (configuration-model {1 #{2 3}
                        2 #{1 5}
                        3 #{1 2 4}
                        4 #{3}
                        5 #{2}}
                       10))


(defn pat-counts [g]
  (let [subgraphs (enumerate-subgraphs g 4)
        {A [3 3 3 3]
         B [2 2 3 3]
         C [2 2 2 2]
         D [1 1 2 2]} (degree-seq-counts subgraphs)]
    [A B C D]))


(def ers-fut
  (future
    (def ers-counts (->> (for [i (range 25)]
                           (do
                             (prn "er" i)
                             (time (pat-counts (erdos-renyi 128 2075)))))
                         doall))))


(def bars-fut
  (future
    (def bars-counts (->> (for [i (range 25)]
                            (do
                              (prn "bars" i)
                              (time (pat-counts (barabasi-albert (erdos-renyi 18 95) 110 18)))))
                          doall))))


(def config-fut
  (future
    (def config-counts (->> (for [i (range 25)]
                              (do
                                (prn "config" i)
                                (time (pat-counts (configuration-model g 1000)))))
                            doall))))


(defn z-score [c cm sd]
  (/ (- c cm)
     sd))


(defn -main [& args]
  @ers-fut @bars-fut @config-fut ;; HACK wait for all futures to finish
  (let [{A [3 3 3 3]
         B [2 2 3 3]
         C [2 2 2 2]
         D [1 1 2 2]} g-degree-counts
        g-count [A B C D]
        ers-mean (s/mean ers-counts)
        ers-sd (s/sd ers-counts)
        bars-mean (s/mean bars-counts)
        bars-sd (s/sd bars-counts)
        config-mean (s/mean config-counts)
        config-sd (s/sd config-counts)]
    (prn
     [[:erdos ers-mean ers-sd (mapv z-score g-count ers-mean ers-sd)]
      [:barabasi bars-mean bars-sd (mapv z-score g-count bars-mean bars-sd)]
      [:config config-mean config-sd (mapv z-score g-count config-mean config-sd)]])))


(comment
  (s/mean ers-counts)

  (s/sd [[1 2 3 4]
         [1 2 3 4]])

  )




