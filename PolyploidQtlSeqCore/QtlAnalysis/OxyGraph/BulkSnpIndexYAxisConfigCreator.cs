namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Bulk SNP-indexグラフのY軸設クリエーター
    /// </summary>
    internal static class BulkSnpIndexYAxisConfigCreator
    {
        private static readonly double _padding = 0.02;
        private static readonly double _step = 0.1;

        /// <summary>
        /// Bulk SNP-indexグラフのY軸設定を作成する。
        /// </summary>
        /// <param name="variants">変異</param>
        /// <returns>Y軸設定</returns>
        public static YAxisConfig Create(SnpIndexVariantWithSlidingWindowQtl[] variants)
        {
            var maxSnpIndex = variants.Select(x => x.Bulk1.SnpIndex.Value)
                .Concat(variants.Select(x => x.Bulk2.SnpIndex.Value))
                .Max();

            var min = -1 * _padding;
            var max = maxSnpIndex + _padding;

            return new YAxisConfig(min, max, _step);
        }
    }
}
