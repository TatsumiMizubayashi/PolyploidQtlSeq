namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// 解析対象変異ルール
    /// </summary>
    internal interface IAnalyzableVariantRule
    {
        /// <summary>
        /// 変異が解析対象かどうかを判断する。
        /// </summary>
        /// <param name="variant">変異</param>
        /// <returns>解析対象ならtrue</returns>
        bool Ok(VcfVariant variant);
    }
}
