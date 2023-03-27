namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// 親2 Depthルール
    /// </summary>
    internal class Parent2DepthRule : IQtlSeqTargetVariantRule
    {
        private readonly MinimumDepthThreshold _threshold;

        /// <summary>
        /// 親2 Depthルールを作成する。
        /// </summary>
        /// <param name="minimumDepthThreshold">しきい値</param>
        public Parent2DepthRule(MinimumDepthThreshold minimumDepthThreshold)
        {
            _threshold = minimumDepthThreshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Parent2.Depth >= _threshold.Value;
        }
    }
}
