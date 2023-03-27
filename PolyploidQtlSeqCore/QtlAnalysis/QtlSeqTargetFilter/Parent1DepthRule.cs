namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// 親1Depthルール
    /// </summary>
    internal class Parent1DepthRule : IQtlSeqTargetVariantRule
    {
        private readonly MinimumDepthThreshold _threhold;

        /// <summary>
        /// 親1Depathルールを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Parent1DepthRule(MinimumDepthThreshold threshold)
        {
            _threhold = threshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Parent1.Depth >= _threhold.Value;
        }
    }
}
