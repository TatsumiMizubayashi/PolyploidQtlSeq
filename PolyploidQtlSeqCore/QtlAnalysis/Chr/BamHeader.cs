namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// BAMファイルのヘッダー情報
    /// </summary>
    internal class BamHeader
    {
        private const string SEQ_TAG = "@SQ";
        private const int NAME_INDEX = 2;
        private const int LENGTH_INDEX = 4;
        private static readonly char[] _spliter = ['\t', ':'];

        private readonly string[] _lines;

        /// <summary>
        /// BAMヘッダーを作成する。
        /// </summary>
        /// <param name="headerLines">ヘッダー行</param>
        public BamHeader(string[] headerLines)
        {
            _lines = headerLines;
        }

        /// <summary>
        /// 染色体コレクションに変換する。
        /// </summary>
        /// <returns>染色体コレクション</returns>
        public Chromosomes ToChromosomes()
        {
            var chrList = new List<Chromosome>();

            var seqLines = _lines.Where(x => x.StartsWith(SEQ_TAG));
            foreach (var line in seqLines)
            {
                var items = line.Split(_spliter);

                var chrName = items[NAME_INDEX];
                var length = int.Parse(items[LENGTH_INDEX]);

                chrList.Add(new Chromosome(chrName, length));
            }

            return new Chromosomes(chrList);
        }
    }
}
