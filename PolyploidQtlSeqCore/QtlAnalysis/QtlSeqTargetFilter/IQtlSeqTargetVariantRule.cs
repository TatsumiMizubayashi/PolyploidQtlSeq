namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// QTL-seq対象変異ルールインターフェイス
    /// </summary>
    internal interface IQtlSeqTargetVariantRule
    {
        /// <summary>
        /// QTL-seq解析対象変異かどうかを判断する。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <returns>解析対象変異の場合はtrue</returns>
        bool Ok(SnpIndexVariant variant);
    }
}
