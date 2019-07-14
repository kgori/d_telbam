import std.stdio;
import bio.std.hts.bam.reader;
import bio.std.hts.bam.read;
import bio.std.hts.bam.referenceinfo;
import bio.std.hts.sam.header;
import bio.std.hts.bam.writer;
import std.algorithm : canFind;
import std.algorithm : sort;
import std.array;
import std.parallelism: TaskPool;
import std.container.rbtree;

auto new_bamwriter_from_template(Reader)(string filename, ref Reader reader, ref TaskPool taskpool, bool create_index=false) {
    // Remember to call scope(exit) writer.finish() on the returned writer
    auto writer = new BamWriter(filename, -1, taskpool);
    if (!create_index) {
        writer.disableAutoIndexCreation();
    }
    writer.writeSamHeader(reader.header);
    writer.writeReferenceSequenceInfo(reader.reference_sequences);
    return writer;
}

bool read_is_less(BamRead a, BamRead b) {
    int a_ref_id = a.ref_id == -1 ? int.max : a.ref_id;
    int b_ref_id = b.ref_id == -1 ? int.max : b.ref_id;
    if (a_ref_id < b_ref_id) return true;
    if (a_ref_id == b_ref_id) {
        if (a.position < b.position) {
            return true;
        }
    }
    return false;
}

bool read_is_same(BamRead a, BamRead b) {
    int a_ref_id = a.ref_id == -1 ? int.max : a.ref_id;
    int b_ref_id = b.ref_id == -1 ? int.max : b.ref_id;
    return (a.ref_id == b.ref_id && a.position == b.position) ? true : false;
}

bool read_is_more(BamRead a, BamRead b) {
    int a_ref_id = a.ref_id == -1 ? int.max : a.ref_id;
    int b_ref_id = b.ref_id == -1 ? int.max : b.ref_id;
    if (a_ref_id > b_ref_id) return true;
    if (a_ref_id == b_ref_id) {
        if (a.position > b.position) {
            return true;
        }
    }
    return false;
}

struct ReadPair (R = BamRead) {
    R read1;
    R read2;

    R maxRead() {
        return read_is_more(this.read1, this.read2) ? this.read1 : this.read2;
    }

    R minRead() {
        return (read_is_less(this.read1, this.read2) || read_is_same(this.read1, this.read2)) ? this.read1 : this.read2;
    }

    bool opEquals(ReadPair o) {
        R this_max = this.maxRead();
        R that_max = o.maxRead();
        return (this_max.ref_id == that_max.ref_id) && (this_max.position == that_max.position);
    }

    int opCmp(ReadPair o) {
        R this_max = this.maxRead();
        R that_max = o.maxRead();
        if (read_is_less(this_max, that_max)) {
            return -1;
        }
        else if (read_is_more(this_max, that_max)) {
            return 1;
        }
        else {
            return 0;
        }
    }
}

void main(string[] args)
{
    auto pool = new TaskPool(2);
    scope(exit) pool.finish();
    auto bam = new BamReader(args[1], pool);
    int total_reads = 0;
    int telbams = 0;
    auto rbt = redBlackTree!string();
    foreach (read; bam.reads()) {
        total_reads += 1;
        if (total_reads % 25000 == 0) {
            writefln("Processed %d reads", total_reads); //Found 29900758 reads; 9089 telreads
        }
        if (canFind(read.sequence, "TTAGGGTTAGGG") || canFind(read.sequence, "CCCTAACCCTAA") && !read.is_supplementary && read.is_paired) {
            rbt.insert(read.name);
            telbams++;
        }
    }
	writefln("Found %d reads; %d telreads", total_reads, telbams);

    auto writer = new_bamwriter_from_template(args[2], bam, pool);
    scope(exit) writer.finish();
    ReadPair!BamRead[string] store;
    writeln("Pairing and sorting reads");


    foreach (read; bam.reads()) {
        if (read.name in rbt && !read.is_supplementary && read.is_paired) {
            //writeln(read.name);
            if (read.name in store) {
                store[read.name].read2 = read;
            } else {
                store[read.name] = ReadPair!BamRead();
                store[read.name].read1 = read;
            }
            //writer.writeRecord(read);
        }
    }
    auto sorted_readpairs = sort!"a < b"(store.values).array;
    foreach (readpair; sorted_readpairs) {
        writer.writeRecord(readpair.maxRead());
        writer.writeRecord(readpair.minRead());
    }
}
