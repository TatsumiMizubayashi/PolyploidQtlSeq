namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 染色体コレクション
    /// </summary>
    internal class Chromosomes
    {
        private readonly IReadOnlyDictionary<string, Chromosome> _chrNameDictionary;

        /// <summary>
        /// 染色体コレクションを作成する。
        /// </summary>
        /// <param name="chromosomes">染色体</param>
        public Chromosomes(IEnumerable<Chromosome> chromosomes)
        {
            _chrNameDictionary = chromosomes.ToDictionary(x => x.Name);
        }

        /// <summary>
        /// 解析対象染色体情報を取得する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>染色体配列</returns>
        public Chromosome[] Get(AnalysisChrOption option)
        {
            if (option.AnalysisChrNames.HasNames) return Get(option.AnalysisChrNames);

            return Get(option.ChrSizeThreshold);
        }

        /// <summary>
        /// 指定した名前と一致する染色体情報を取得する。
        /// </summary>
        /// <param name="analysisChr">染色多名</param>
        /// <returns>染色体配列</returns>
        public Chromosome[] Get(AnalysisChrNames analysisChr)
        {
            var chrList = new List<Chromosome>();

            foreach (var name in analysisChr.Names)
            {
                if (!_chrNameDictionary.TryGetValue(name, out var chr))
                    throw new ArgumentException($"Chromosome [{name}] do not exist.");

                chrList.Add(chr);
            }

            return chrList.ToArray();
        }

        /// <summary>
        /// しきい値以上の長さを持つ染色体情報を取得する。
        /// </summary>
        /// <param name="sizeThreshold">長さのしきい値</param>
        /// <returns>染色体配列</returns>
        public Chromosome[] Get(ChrSizeThreshold sizeThreshold)
        {
            return _chrNameDictionary.Values
                .Where(x => x.Length >= sizeThreshold.Value)
                .ToArray();
        }
    }
}
