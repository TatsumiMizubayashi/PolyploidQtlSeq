namespace Sequence.Position.Extract
{
    /// <summary>
    /// GenomePositionを持つ項目のインターフェイス
    /// </summary>
    public interface IHasGenomePositionItem
    {
        /// <summary>
        /// ゲノム位置情報を取得する。
        /// </summary>
        GenomePosition GenomePosition { get; }
    }
}
