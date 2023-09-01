namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Mapping設定値インターフェース
    /// </summary>
    public interface IMappingSettingValue
    {
        /// <summary>
        /// リファレンスシークエンスファイルを取得する。
        /// </summary>
        string ReferenceSequence { get; }

        /// <summary>
        /// スレッド数を取得する。
        /// </summary>
        int ThreadNumber { get; }
    }
}
