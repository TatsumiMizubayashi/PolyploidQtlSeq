namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// Bulk2 最大SNP-indexルール
    /// </summary>
    internal class Bulk2MaxSnpIndexRule : IQtlSeqTargetVariantRule
    {
        private readonly MaxBulkSnpIndexThreshold _threshold;

        /// <summary>
        /// Bulk2 最大SNP-indexルール インスタンスを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Bulk2MaxSnpIndexRule(MaxBulkSnpIndexThreshold threshold)
        {
            _threshold = threshold;
        }

        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Bulk2.SnpIndex.Value <= _threshold.Value;
        }
    }
}
