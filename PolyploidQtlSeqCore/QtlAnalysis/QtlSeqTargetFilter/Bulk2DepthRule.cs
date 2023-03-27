namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Bulk2 Depthルール
    /// </summary>
    internal class Bulk2DepthRule : IQtlSeqTargetVariantRule
    {
        private readonly MinimumDepthThreshold _threshold;

        /// <summary>
        /// Bulk2 Depthルールを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Bulk2DepthRule(MinimumDepthThreshold threshold)
        {
            _threshold = threshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Bulk2.Depth >= _threshold.Value;
        }
    }
}
