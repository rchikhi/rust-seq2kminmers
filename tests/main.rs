fn main() {

    let genome_file = "../test/ecoli.genome.100k.fa";
    let contents = fs::read_to_string(genome_file).expect("Failed to read test genome file").split("\n").collect()[1].strip();

    let kminmers = [331915013485775537, 145328730389169902, 535416482282343629, 310285671002660864];
    let pos = [99843, 99846, 99860, 99884];

    println!("");
}
