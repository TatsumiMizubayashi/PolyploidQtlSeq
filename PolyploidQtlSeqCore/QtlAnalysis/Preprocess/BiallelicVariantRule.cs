namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// bialleleic変異ルール
    /// </summary>
    internal class BiallelicVariantRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return !variant.IsMultiAltAllele;
        }
    }
}
