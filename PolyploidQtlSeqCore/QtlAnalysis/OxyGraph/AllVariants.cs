namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// 全変異
    /// </summary>
    internal class AllVariants
    {
        private readonly string[] _chrNames;
        private readonly SnpIndexVariantWithSlidingWindowQtl[] _variants;
        private readonly IReadOnlyDictionary<string, ChrVariants> _chrVariantsDictionary;

        /// <summary>
        /// 全変異コレクションを作成する。
        /// </summary>
        /// <param name="variants">全変異</param>
        public AllVariants(SnpIndexVariantWithSlidingWindowQtl[] variants)
        {
            _variants = variants;
            _chrNames = _variants.Select(x => x.GenomePosition.ChrName).Distinct().ToArray();
            _chrVariantsDictionary = _variants.ToLookup(x => x.GenomePosition.ChrName)
                .ToDictionary(x => x.Key, x => new ChrVariants(x.Key, x.ToArray()));
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        /// <returns>染色体名</returns>
        public string[] GetChrNames() => _chrNames;

        /// <summary>
        /// 指定染色体の変異コレクションを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>指定染色体の変異コレクション</returns>
        public ChrVariants GetVariants(string chrName) => _chrVariantsDictionary[chrName];

        /// <summary>
        /// Bulk SnpIndexのY軸設定を作成する。
        /// </summary>
        /// <returns>Bulk SnpIndexのY軸設定</returns>
        public YAxisConfig CreateBulkSnpIndexYAxisConfig()
        {
            return BulkSnpIndexYAxisConfigCreator.Create(_variants);
        }
    }
}
