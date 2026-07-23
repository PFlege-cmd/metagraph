//
// Created by Patrick Flege on 14/07/2026.
//

#include "RegionOfInterestAdapter.h"

#include "CoordinateRetriever.h"
#include "GenomeCoordinateFactory.h"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "graph/annotated_dbg.hpp"

struct SequencePositions {
    int no_of_results;
    int* sequences;
    int* positions;
    SequencePositions() = default;
    SequencePositions(int no_results, int* seqs, int* poss ) { no_of_results = no_results; sequences = seqs; positions = poss; };
};
extern "C" {
    SequencePositions* retrieveKmersOfInterest(char *database, int genome_nr, char * genome_name, int no_of_kmers, int no_of_sequences, const char ** kmers, int * sequences) {

        graph_glue glue = graph_glue(0, NULL);
        std::vector<int> lengths = std::vector<int>(no_of_sequences );
        std::vector<std::string_view> kmer_views = std::vector<std::string_view>( no_of_kmers);
        std::vector<int> idxs = std::vector<int>(no_of_kmers);
        std::vector<int> seqs_idxs = std::vector<int>(no_of_sequences);

        std::iota(seqs_idxs.begin(), seqs_idxs.end(), 0);
        for_each(seqs_idxs.begin(), seqs_idxs.end(), [&sequences, &lengths](int idx) {
            lengths[idx] = sequences[idx];
        });

        std::iota(idxs.begin(), idxs.end(), 0);
        for_each(idxs.begin(), idxs.end(), [&kmer_views, &kmers](int idx) {
            auto kv = std::string_view(kmers[idx]);
            kmer_views[idx] = std::string_view(kv);
        });

        std::vector<unique_ptr<MatrixEntry>> coordinates
            = std::vector<std::unique_ptr<MatrixEntry>>();

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(database);
        auto wrapper = DeBruijnGraphWrapper(*graph);

        static std::shared_ptr<AnnotatedDBG> kmer_graph = glue.load_dbg(database);
        auto kmer_graph_wrapper = DeBruijnGraphWrapper(*kmer_graph);

        auto kmer_labels = kmer_graph_wrapper.get_graph()->get_annotator().get_label_encoder().get_labels();
        auto coord_labels = wrapper.get_graph()->get_annotator().get_label_encoder().get_labels();

        std::string genome_number_string = kmer_labels[genome_nr];
        std::string_view genome_name_view = genome_name;

        std::cout << "Kmer label is: " <<  genome_number_string << std::endl;
        std::cout << "Coordinate label is: " <<  genome_name_view << std::endl;

        CoordinateRetriever* retriever = new CoordinateRetriever(no_of_sequences, 5, genome_number_string, sequences, kmer_graph_wrapper);
        kmer_views = retriever->filter_max_frequency(kmer_views);

        for_each(kmer_views.begin(), kmer_views.end(), [&coordinates, &genome_name_view, &wrapper, &lengths](auto& kmer) {
            coordinates.push_back(std::make_unique<GenomeCoordinate>(genome_name_view, kmer, lengths, wrapper));
        });

        auto pos = retriever->get_kmer_positions(coordinates);

        int * found_sequences = new int[pos.size()];
        int * found_positions = new int[pos.size()];
        int number_results = (int) pos.size();

        for (int i = 0; i < number_results; i++) {
            found_sequences[i] = pos[i][0];
            found_positions[i] = pos[i][1];
        }
        auto final_res = new SequencePositions(number_results, found_sequences, found_positions);
        return final_res;
    }

}