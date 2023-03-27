using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// ΔSNP-index Y軸設定クリエーター
    /// </summary>
    internal static class DeltaSnpIndexYAxisConfigCreator
    {
        private static readonly double _padding = 0.02;
        private static readonly double _step = 0.1;

        /// <summary>
        /// ΔSNP-index Y軸設定を作成する。
        /// </summary>
        /// <param name="variants">変異</param>
        /// <param name="windows">window</param>
        /// <returns>ΔSNP-index Y軸設定</returns>
        public static YAxisConfig Create(SnpIndexVariantWithSlidingWindowQtl[] variants, Window[] windows)
        {
            var minDeltaSnpIndex = variants.Min(x => x.DeltaSnpIndex.Value);
            var maxDeltaSnpIndex = variants.Max(x => x.DeltaSnpIndex.Value);
            var maxP99Threshold = windows.Max(x => x.P99Qtl.ThresholdDeltaSnpIndex.Value);

            var values = new[] { minDeltaSnpIndex, maxDeltaSnpIndex, maxP99Threshold, -1 * maxP99Threshold };
            var min = values.Min() - _padding;
            var max = values.Max() + _padding;

            return new YAxisConfig(min, max, _step);
        }
    }
}
