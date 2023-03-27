namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Bulk1 Depthルール
    /// </summary>
    internal class Bulk1DepthRule : IQtlSeqTargetVariantRule
    {
        private readonly MinimumDepthThreshold _threshold;

        /// <summary>
        /// Bulk1 Depthルールを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Bulk1DepthRule(MinimumDepthThreshold threshold)
        {
            _threshold = threshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Bulk1.Depth >= _threshold.Value;
        }
    }
}
